
#include <cuda_runtime.h>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <fcntl.h>
#include <unistd.h>
#include <sched.h>
#include <pthread.h>
#include "sem.h"
#include <errno.h>
#include <time.h>
#include "sampleblock.h"

#define SAMPLEBLOCK_THREAD_PRIORITY     28

unsigned char dsp::SampleBlock::NumInstances;

pthread_mutex_t dsp::SampleBlock::mtx;

dsp::SampleBlock::SampleBlock(){
    if (NumInstances == 0)
        pthread_mutex_init(&mtx, NULL);
    pthread_mutex_lock(&mtx);
    NumInstances++;
    pthread_mutex_unlock(&mtx);
    fd = -1;
    ModuleName = "SampleBlock";
    AllocateInputs(1);
    AllocateOutputs(3);
    KeepRunning = false;
    InputSourceType = SAMPLE_INPUT_SOURCE_FILE;
    
    
    ConfigExpectedInput(0, "StartByte", INT_t, VALUE, 1);

    // VectorLength to be initialised at Open() Set to 2 to prevent errors for now
    // *Data will be updated at every call to Update()
    ConfigOutput(0, "Samples", UNDEFINED_t, VALUE_CMPX, CUDA_DEVICE, 2, NULL, 0);

    ConfigOutput(1, "SamplingFrequency", DOUBLE_t, FREQUENCY_HZ, HOST, 1, (void*)&SamplingFrequency, 0);
    ConfigOutput(2, "SampleLength", DOUBLE_t, VALUE, HOST, 1, (void*)&SampleLength, 0);

    std::clog << "[" << ModuleName << "] Configured outputs" << std::endl;

    cudaStreamCreate(&cuStream);
    
    InsertParam("Filename", (void*)&Filename, CHAR_t, FilenameCapacity, 0);
    InsertParam("Hostname", (void*)&Hostname, CHAR_t, HostnameCapacity, 0); // enter as "192.168.10.7"
    InsertParam("PortNo", (void*)&PortNo, INT_t, sizeof(int), sizeof(int)); // enter as 49152
    InsertParam("SamplingFrequency", (void*)&SamplingFrequency, DOUBLE_t,
                sizeof(double), sizeof(double)); // This is in Hz
    InsertParam("SampleLength", (void*)&SampleLength, DOUBLE_t, sizeof(double), sizeof(double));
    InsertParam("RunLive", (void*)&RunLive, BOOL_t, sizeof(bool), sizeof(bool));
    InsertParam("InputSourceType", (void*)&InputSourceType, CHAR_t,
                sizeof(char), sizeof(char)); // Probably need to change this... need to enter as \x1 or \x0
}

dsp::SampleBlock::~SampleBlock(){
    if (KeepRunning) Stop();
    //Close();
    pthread_mutex_lock(&mtx);
    NumInstances--;
    pthread_mutex_unlock(&mtx);
    delete [] outputs;
    cudaStreamDestroy(cuStream);
}

int dsp::SampleBlock::Close() {
    if (fd < 0) {
        std::clog << "[" << ModuleName << "] No file is open." << std::endl;
        return -1;
    }
    if (KeepRunning) {
        std::cerr << "[" << ModuleName << "] Can't close file while running."
                  << std::endl;
        return -1;
    }
    int i = close(fd); fd = 0;
    if (i) {
        std::cerr << "[" << ModuleName << "] Error closing file.";
        return -1;
    }
    for (int i = 0; i < NumBlocks; i++) {
        cudaFreeHost(BlockArr[i]);
        cudaFree(BlockArrCuda[i]);
    }
    std::free(BlockArrCuda);
    std::free(BlockArr);
    NumBlocks = 0;
    SamplingFrequency = 0.0;
    this->BlockLength = 0;
    CurrLoadBlock = -1;
    CurrProcBlock = -1;
    outputs[0].VectorLength = 0;
    outputs[0].Data = NULL;
    return 0;
}

int dsp::SampleBlock::Open() {
    int i, j;

    if (fd >= 1){
        std::clog << "[SampleBlock] An existing file is already open."
                  << std::endl;
        return -1;
    }

    // TODO: consider cudaMallocHost for memory that needs copied between device and host regularly (samples)

    switch (InputSourceType){
    case SAMPLE_INPUT_SOURCE_FILE:
        // Open the file (using C operations for speed?)
        fd = open(this->Filename, O_RDONLY);
        if (fd < 0) {
            std::cerr << "[SampleBlock] Unable to open file: " << this->Filename
                      << std::endl;
            return -1;
        }
        // Skip ahead the number of bits read by the scalar tracker
        skipCheck = lseek(fd, *((long long int*)(inputs[0]->Data)), SEEK_SET);
        if (skipCheck == 0) {
            std::cerr << "[SampleBlock] Failed to skip ahead in file: " << this->Filename
                      << std::endl;
            return -1;
        }
        else {
        	std::clog << "[" << ModuleName << "] Starting reading at byte " << lseek(fd, 0, SEEK_CUR)
        		      << " in file " << this->Filename << std::endl;
        }
        break;
    case SAMPLE_INPUT_SOURCE_SOCKET:
        struct hostent *server;
        memset((void*)&ServerAddr, 0, sizeof(sockaddr_in));
        fd = socket(AF_INET, SOCK_STREAM, 0);
        if (fd < 0) {
            std::cerr << "[" << ModuleName << "] Open: Unable to open socket." << std::endl;
            return -1;
        }
        server = gethostbyname(Hostname);
        if (server == NULL){
            std::cerr << "[" << ModuleName << "] Open: Invalid hostname." << std::endl;
            close(fd); fd = 0;
            return -1;
        }
        ServerAddr.sin_family = AF_INET;
        memcpy((void*)&ServerAddr.sin_addr.s_addr, (void*)server->h_addr, server->h_length);
        ServerAddr.sin_port = htons(PortNo);
        if (connect(fd,(struct sockaddr *)&ServerAddr,sizeof(ServerAddr)) < 0){
            std::cerr << "[" << ModuleName << "] Open: Failed to connect." << std::endl;
            close(fd); fd = 0;
            return -1;
        }
        break;
    default:
        std::cerr << "[" << ModuleName << "] Open: Invalid input source type."
                  << std::endl;
        break;
    }

    NumBlocks = NumBlocksDefault;
    // Sets the file read to be 1ms' worth
    // (fs/1000 = # of samples in 1ms = samples per mHz)
    // (The +0.5 rounds, since BlockLength is an int, so if 1ms is more than half way
    //  to the next sample, include it)
    //this->BlockLength = (unsigned short)(SamplingFrequency / 1000.0 + 0.5);
    this->BlockLength = (unsigned short)(SamplingFrequency * SampleLength + 0.5); // Attempt at custom-sized sample lengths
    CurrLoadBlock = -1;
    CurrProcBlock = -1;

    BlockArr = reinterpret_cast<void **>(std::malloc(sizeof(void *) * NumBlocks));
    if (BlockArr == NULL) {
        close(fd); fd = 0;
        NumBlocks = 0;
        this->BlockLength = 0;
        std::cerr << "[SampleBlock] Unable to allocate sample buffers.\n"
                  << "File: " << __FILE__ << " Line: " << __LINE__
                  << std::endl;
        return -1;
    }

    BlockArrCuda = reinterpret_cast<void **>(std::malloc(sizeof(void *) * NumBlocks));
    if (BlockArr == NULL) {
        std::free(BlockArr);
        close(fd); fd = 0;
        NumBlocks = 0;
        this->BlockLength = 0;
        std::cerr << "[SampleBlock] Unable to allocate cuda device sample buffers.\n"
                  << "File: " << __FILE__ << " Line: " << __LINE__
                  << std::endl;
        return -1;
    }

    size_t blockSizeBytes = sizeof(int16_t)*BlockLength*2;

    for (i=0; i<NumBlocks; i++) {
        /* Memory aligned to 64 for cache speedups and future SIMD uses.
           Aligned to 4096 to get its own page in virtual memory -->
           TLB speedups tho negligible */
         //j = posix_memalign(&BlockArr[i], 4096, sizeof(int16_t)*BlockLength*2);
         //if (j < 0) {
         cudaError_t err;
         err = cudaMallocHost(&BlockArr[i], blockSizeBytes);
         if (err != cudaSuccess){
            std::cerr << "[SampleBlock] Unable to allocate sample buffers.\n"
                      << "File: " << __FILE__ << " Line: " << __LINE__
                      << "\ni = " << i << "\n" << cudaGetErrorString(err)
                      << std::endl;
            for (j = 0; j < i; j++) {
                cudaFreeHost(BlockArr[j]);
                cudaFree(BlockArrCuda[j]);
            }
            std::free(BlockArr);
            close(fd); fd = 0;
            NumBlocks = 0;
            this->BlockLength = 0;
            return -1;
         }
         err = cudaMalloc(&BlockArrCuda[i], blockSizeBytes);
         if (err != cudaSuccess) {
            std::cerr << "[SampleBlock] Unable to allocate cuda device sample buffers.\n"
                      << "File: " << __FILE__ << " Line: " << __LINE__
                      << "\ni = " << i << "\n" << cudaGetErrorString(err)
                      << std::endl;
            std::free(BlockArr[i]);
            for (j = 0; j < i; j++) {
                cudaFreeHost(BlockArr[j]);
                cudaFree(BlockArrCuda[j]);
            }
            std::free(BlockArrCuda);
            std::free(BlockArr);
            close(fd); fd = 0;
            NumBlocks = 0;
            this->BlockLength = 0;
            return -1;
         }
    }

    // Note: the vector length is the number of complex samples,
    //       meaning the actual data size is 2*BlockLength*sizeof(type)
    //this->RealTime = RunLive;
    outputs[0].VectorLength = BlockLength;

    return 0;
}

int dsp::SampleBlock::Start(void* cuFlowStream){
    int i;
    pthread_attr_t attr;
    struct sched_param param;
    /*if (fd < 0) {
        std::cerr << "[SampleBlock] File not open." << std::endl;
        return -1;
    }*/
    std::clog << "[" << ModuleName << "] Starting ... " << std::flush;
    i = Open();
    if (i) {
        std::cerr << "[SampleBlock] File not open." << std::endl;
        return -1;
    }
    i = pthread_attr_init(&attr);
    if (i) {
        std::cerr << "[SampleBlock] Unable to init thread attr."
                  << std::endl;
        return -1;
    }
    int minPrio = sched_get_priority_min(SCHED_RR);
    int maxPrio = sched_get_priority_max(SCHED_RR);
    /* Posix min requirements is 32 priority levels, Although linux provides
       99 levels, good to check for min/max available levels. */
    param.sched_priority = (maxPrio - minPrio) *
                           SAMPLEBLOCK_THREAD_PRIORITY / 32 + minPrio;
    i = pthread_attr_setschedpolicy(&attr, SCHED_RR);
    if (i) {
        std::cerr << "[SampleBlock] Unable to set sample thread real-time."
                  << std::endl;
        return -1;
    }
    i = pthread_attr_setschedparam(&attr, &param);
    if (i) {
        std::cerr << "[SampleBlock] Unable to set sample thread priority."
                  << std::endl;
        return -1;
    }
    i = pthread_create(&thread, &attr, SampleThreadEntry, this);
    if (i) {
        std::cerr << "[SampleBlock] Unable to create thread."
                  << std::endl;
        return -1;
    }
    while (KeepRunning == false);
    std::clog << "Started." << std::endl;
    return 0;
}

int dsp::SampleBlock::Stop(){
    KeepRunning = false;
    FirstUpdate = false;
    std::clog << "[" << ModuleName << "] Stopping ... " << std::flush;
    int ret = pthread_join(thread, NULL);
    std::clog << "Stopped." << std::endl;
    return ret;
}

void * dsp::SampleBlock::SampleThreadEntry(void * thisObj){
    ((dsp::SampleBlock *)thisObj)->GetSamplesThread();
    return NULL;
}

void dsp::SampleBlock::GetSamplesThread(){
	// Set up semaphores and buffer tracking indices
    ssize_t BytesRead;
    cudaEvent_t cuEvent;
    uint8_t firstIteration = 1;
    int8_t reachedEOF = 0;
    uint64_t blockCnt = 0;
    if (fd < 0) {
        std::cerr << "[SampleBlock] File not open. Quitting thread."
                  << std::endl;
        return;
    }
    CurrLoadBlock = 0;
    CurrLoadCudaBlock = 0;
    CurrProcBlock = -1;
    if (osindep_sem_init(&samplesAvailSem, 0, 0)){
        std::cerr << "[SampleBlock] Can't initialize semaphores. Quiting thread."
                  << std::endl;
        return;
    }
    if (osindep_sem_init(&buffersAvailSem, 0, (NumBlocks-1))){
        osindep_sem_destroy(&samplesAvailSem);
        std::cerr << "[SampleBlock] Can't initialize semaphores. Quiting thread."
                  << std::endl;
        return;
    }
    cudaEventCreate(&cuEvent);
    FirstUpdate = true;
    KeepRunning = true;

    // Sit and spin until the first update happens so:
    //  - The semaphore timeouts don't expire in post-processing mode
    // or
    //  - The newest samples are fetched once the pipeline is ready to run in RunLive mode
    // Note: Update will block until one buffer holds data, so it's safe to not preload data until FirstUpdate is underway
    while (FirstUpdate)
        usleep(200);
    while (KeepRunning){
        if (!reachedEOF){
            BytesRead = 0;
            int readCnt;
            // Read until the entire 1ms block has been read in
            // ("read" deals in bytes read, so 4*BlockLength accounts for
            //  a real and imaginary sample, each int16_t = 2bytes read in)
            while (BytesRead < (4*BlockLength)){
                void *ptr = (void*)&((char*)BlockArr[CurrLoadBlock])[BytesRead];
                readCnt = read(fd, ptr, (4*BlockLength-BytesRead)); // Note this is C's "read" function... why??? speed?
                if (readCnt < 0){
                    KeepRunning = false;
                    std::cerr << "[SampleBlock] Read Error. File: "
                              << __FILE__ << " Line: " << __LINE__
                              << " CurrLoadBlock: " << CurrLoadBlock
                              << std::endl;
                    perror(NULL);
                    break;
                } else if (readCnt == 0){
                    KeepRunning = false;
                    reachedEOF = 1;
                    std::clog << "[SampleBlock] Reached EOF." << std::endl;
                    std::clog << "[SampleBlock] blockCnt = " << blockCnt << std::endl;
                    break;
                }
                BytesRead += readCnt;
            }
        }

        if (firstIteration){
            firstIteration = 0;
        } else {
            if (cudaEventSynchronize(cuEvent) == cudaSuccess){
                CurrLoadCudaBlock++;
                if (CurrLoadCudaBlock >= NumBlocks) CurrLoadCudaBlock = 0;

                if (osindep_sem_post(&samplesAvailSem)){
                    std::cerr << "[SampleBlock] Error: sem_post in file: "
                              << __FILE__ << " line: " << __LINE__
                              << std::endl;
                    KeepRunning = false;
                    return;
                }
            } else {
                KeepRunning = false;
                std::cerr << "[SampleBlock] Error: cudaEventSynchronize in file: "
                          << __FILE__ << " line: " << __LINE__
                          << std::endl;
            }
        }

        if (BytesRead == (4*BlockLength)){
        //    if (cudaMemcpy(BlockArrCuda[CurrLoadBlock], BlockArr[CurrLoadBlock],
        //            BytesRead, cudaMemcpyHostToDevice) != cudaSuccess){
            if (cudaMemcpyAsync(BlockArrCuda[CurrLoadCudaBlock], BlockArr[CurrLoadBlock],
                    BytesRead, cudaMemcpyHostToDevice, cuStream) != cudaSuccess){
                KeepRunning = false;
                std::cerr << "[SampleBlock] Error: cudaMemcpy in file: "
                          << __FILE__ << " line: " << __LINE__ << std::endl;
            } else {
                cudaEventRecord(cuEvent, cuStream);
            }
            blockCnt++;
        }

        if (reachedEOF) break;

        CurrLoadBlock++;
        if (CurrLoadBlock >= NumBlocks) CurrLoadBlock = 0;


        // Check if there are any buffers available (without blocking)
        if (osindep_sem_trywait(&buffersAvailSem)){

        	// In RunLive mode, if no buffers are available immediately, samples are going to be lost, so crash
            if (RunLive)
                std::clog << "[SampleBlock] Fail real-time." << std::endl;


            // In post-process mode, wait for a reasonable amount of time for the buffers to free up
            // before assuming something got stuck somewhere.

            // timeout after 1.5s
            if (osindep_sem_waitforduration(&buffersAvailSem,  1500000000)){
			//if (osindep_sem_waitforduration(&buffersAvailSem, 10000000000000)){ // TODO: Remove -- extended timeout for debugging
                if (errno == ETIMEDOUT)
                    std::cerr << "[SampleBlock] Error: sem_timewait timeout: buffersAvailSem"
                              << std::endl;
                else if (errno == EINTR)
                    std::cerr << "[SampleBlock] Error: sem_timewait EINTR: buffersAvailSem"
                              << std::endl;
                else
                    std::cerr << "[SampleBlock] Error: sem_timewait: buffersAvailSem"
                              << std::endl;
                osindep_sem_destroy(&samplesAvailSem);
                osindep_sem_destroy(&buffersAvailSem);
                KeepRunning = false;
                return;
            }
        }
    }

    // Wait till Flow has finished processing all loaded samples
    int semVal = 1;
    while ((semVal) && (KeepRunning)){
        if (osindep_sem_getvalue(&samplesAvailSem, &semVal)){
            std::cerr << "[SampleBlock] Error: sem_getvalue." << std::endl;
            break;
        }
        usleep(1000);
    }
    KeepRunning = false;
    cudaEventDestroy(cuEvent);
    Close();
}

int dsp::SampleBlock::Update(void* cuFlowStream){
    if (KeepRunning) {

    	// Don't sem_post on the first call of Update() --
    	// need one full pipeline iteration before the first block can be safely updated
    	// so let the pipeline get (and stay) one buffer ahead
        if (FirstUpdate) {
        	FirstUpdate = false;
        }
        else {
        	// sem_post so the buffer in the LAST iteration can be replaced
        	if (osindep_sem_post(&buffersAvailSem)){
				std::cerr << "[SampleBlock] sem_post." << std::endl;
				return -1;
			}
        }

        // Block the host until there is a buffer with data available to the pipeline.
        // 1.5ms timeout
        if (osindep_sem_waitforduration(&samplesAvailSem,  1500000000)){
		//if (osindep_sem_waitforduration(&samplesAvailSem, 5000000000000)){ // Remove -- extended timeout for debugging
            if (KeepRunning == 0)
                return -1;
            if (errno == ETIMEDOUT)
                std::cerr << "[SampleBlock] Error: sem_timewait timeout: samplesAvailSem"
                          << std::endl;
            else if (errno == EINTR)
                std::cerr << "[SampleBlock] Error: sem_timewait EINTR: samplesAvailSem"
                          << std::endl;
            else
                std::cerr << "[SampleBlock] Error: sem_timewait: samplesAvailSem"
                          << std::endl;
            osindep_sem_destroy(&samplesAvailSem);
            osindep_sem_destroy(&buffersAvailSem);
            KeepRunning = false;
            return -1;
        }

		// Advance the index of block currently being processed so the output looks at the right buffer
		CurrProcBlock++;
		if (CurrProcBlock >= NumBlocks) CurrProcBlock = 0;

        // update output to point to the right buffer
        outputs[0].Data = BlockArrCuda[CurrProcBlock];

        return 0;
    } else {
        std::cerr << "[SampleBlock] Thread not running." << std::endl;
        return -1;
    }
}


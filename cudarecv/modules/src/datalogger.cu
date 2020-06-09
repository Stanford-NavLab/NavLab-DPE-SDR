
#include <cuda_runtime.h>
#include <cuComplex.h>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>
#if defined(HAS_OPENMP)
    #include <omp.h>
#endif
#include "sem.h"
#include "datalogger.h"
#include "errorhandler.h"

#define cuCheckOpen(stmt)       cuCheckHelper(stmt, Close(); return -1)

#define DATALOGGER_THREAD_PRIORITY     8

dsp::DataLogger::DataLogger(std::string Name) {
    ModuleName = Name;
    ConstructorHelper();
}

dsp::DataLogger::DataLogger() {
    ModuleName = "DataLogger";
    ConstructorHelper();
}

void dsp::DataLogger::ConstructorHelper() {
    AllocateInputs(1);
    AllocateOutputs(0);

    ConfigExpectedInput(0, "Data", DATATYPE_ANY, VALUETYPE_ANY, VECTORLENGTH_ANY);
    InsertParam("Filename", (void*)&Filename, CHAR_t, FilenameCapacity, 0);
    InsertParam("CSV", (void*)&csv, BOOL_t, sizeof(bool), sizeof(bool));
}

dsp::DataLogger::~DataLogger() {
    delete [] expectedInputs;
    delete [] inputs;
    delete [] outputs;
}

int dsp::DataLogger::Open() {
    if (csv)    fd = fopen(Filename, "w");
    else        fd = fopen(Filename, "wb");
    if (fd == NULL) {
        std::cerr << "[" << ModuleName << "] File " << Filename
                  << " cannot be opened." << std::endl;
        return -1;
    }
    cuCheckOpen(cudaMallocHost(&buf_h, bufSize*numBufs));
    for (int i=0; i<numBufs; i++)
        bufEnd[i] = 0;
    currUpdateBuf = 0;
    currThreadBuf = 0;
    return 0;
}

int dsp::DataLogger::Close() {
    if (fd == NULL) {
        std::clog << "[" << ModuleName << "] Nothing to close." << std::endl;
        return 0;
    }
    fclose(fd);
    fd = NULL;
    cuCheck(cudaFreeHost(buf_h));
    return 0;
}

int dsp::DataLogger::Start(void* cuFlowStream) {
    pthread_attr_t attr;
    struct sched_param param;
    std::clog << "[" << ModuleName << "] Starting ... " << std::flush;

    errCheckMod(Open());

    errCheckMod(pthread_attr_init(&attr));
    int minPrio = sched_get_priority_min(SCHED_RR);
    int maxPrio = sched_get_priority_max(SCHED_RR);
    /* Posix min requirements is 32 priority levels, Although linux provides
       99 levels, good to check for min/max available levels. */
    param.sched_priority = (maxPrio - minPrio) *
                           DATALOGGER_THREAD_PRIORITY / 32 + minPrio;
    errCheckMod(pthread_attr_setschedpolicy(&attr, SCHED_RR));
    errCheckMod(pthread_attr_setschedparam(&attr, &param));
    if (pthread_create(&thread, &attr, LoggerThreadEntry, this)) {
        std::cerr << "[" << ModuleName << "] Error: Unable to create thread."
                  << std::endl;
        return -1;
    }
    while (KeepRunning == false);
    std::clog << "Started." << std::endl;
    return 0;
}

int dsp::DataLogger::Stop() {
    std::clog << "[" << ModuleName << "] Stopping ... " << std::flush;
    KeepRunning = false;
    // signal thread to write remaining buffer and quit
    errCheckMod(osindep_sem_post(&bufFullSem));
    int ret = pthread_join(thread, NULL);
    std::clog << "Stopped." << std::endl;
    return ret;
}

void * dsp::DataLogger::LoggerThreadEntry(void * thisObj){
    ((dsp::DataLogger *)thisObj)->LoggerThread();
    return NULL;
}

void dsp::DataLogger::LoggerThread() {
    if (osindep_sem_init(&bufFullSem, 0, 0)) {
        std::cerr << "[" << ModuleName << "] Cannot initialise semaphore." << std::endl;
        return;
    }
    KeepRunning = true;
    while (KeepRunning) {
        // timeout after 1.5s
        if (osindep_sem_waitforduration(&bufFullSem, 1500000000)){
    	//if (osindep_sem_waitforduration(&bufFullSem, 15000000000000)){ // Debug x10000
            if (errno == ETIMEDOUT)
                std::cerr << "[" << ModuleName
                          << "] Error: sem_timewait timeout: buffersAvailSem" << std::endl;
            else if (errno == EINTR)
                std::cerr << "[" << ModuleName
                          << "] Error: sem_timewait EINTR: buffersAvailSem" << std::endl;
            else
                std::cerr << "[" << ModuleName
                          << "] Error: sem_timewait: buffersAvailSem" << std::endl;
            KeepRunning = false;
            break;
        }

        // Get the pointer to the buffer to write (since buffers are only segmented by the way pointers handle things)
        void *ptr = (void*)((char*)buf_h + currThreadBuf*bufSize);

        // Determine how many elements exist in this buffer (due to VectorLengths, the buffer may not be totally full)
        size_t size = bufEnd[currThreadBuf];
        if (csv) {
            switch (inputs[0]->Datatype){
            case FLOAT_t:   size /= sizeof(float);  break;
            case DOUBLE_t:  size /= sizeof(double); break;
            case CHAR_t:    size /= sizeof(char);   break;
            case INT_t:     size /= sizeof(int);    break;
            case BOOL_t:    size /= sizeof(bool);   break;
            default:        continue;   // other types no supported
            }

            // The divide-by-2 isn't needed here since bufEnd keeps track of how full the buffer is in bytes
            //if (inputs[0]->ValueType == VALUE_CMPX)
            //    size >>= 1; //divide by 2

            // Format the output based on datatype
            int col = 0;
            size_t idx = 0;
            while (idx < size) {
                switch (inputs[0]->Datatype){
                case FLOAT_t:   fprintf(fd, "%f", ((float*)ptr)[idx]);  break;
                case DOUBLE_t:  fprintf(fd, "%f", ((double*)ptr)[idx]);  break;
                case CHAR_t:    fprintf(fd, "%d", ((char*)ptr)[idx]);  break;
                case INT_t:     fprintf(fd, "%d", ((int*)ptr)[idx]);  break;
                case BOOL_t:    fprintf(fd, "%d", ((bool*)ptr)[idx]);  break;
                default:        break;   // other types no supported
                }
                if (inputs[0]->ValueType == VALUE_CMPX) {
                    idx++;
                    if (idx >= size) break;
                    switch (inputs[0]->Datatype){
                    case FLOAT_t: {
                        float val = ((float*)ptr)[idx];
                        if (val >= 0.0) fprintf(fd, "+%fj", val);
                        else fprintf(fd, "%fj", val);
                        } break;
                    case DOUBLE_t: {
                        double val = ((double*)ptr)[idx];
                        if (val >= 0.0) fprintf(fd, "+%fj", val);
                        else fprintf(fd, "%fj", val);
                        } break;
                    case CHAR_t: {
                        signed char val = ((char*)ptr)[idx];
                        if (val >= 0) fprintf(fd, "+%dj", val);
                        else fprintf(fd, "%dj", val);
                        } break;
                    case INT_t: {
                        signed int val = ((int*)ptr)[idx];
                        if (val >= 0) fprintf(fd, "+%dj", val);
                        else fprintf(fd, "%dj", val);
                        } break;
                    case BOOL_t:    break;  //bool should not be complex
                    default:        break;  // other types no supported
                    }
                }
                idx++;
                col++;
                if (col >= inputs[0]->VectorLength) {
                    fprintf(fd, "\n");
                    col = 0;
                } else {
                    fprintf(fd, ", ");
                }
            }
        } else {
            fwrite(ptr, 1, size, fd);
        }
        bufEnd[currThreadBuf] = 0;
        currThreadBuf++;
        if (currThreadBuf >= numBufs) currThreadBuf = 0;
    }
    osindep_sem_destroy(&bufFullSem);
    Close();
}

int dsp::DataLogger::Update(void* cuFlowStream) {
    if (KeepRunning) {

    	// Determine how much space to allocate for each input element
        size_t size;
        switch (inputs[0]->Datatype){
        case FLOAT_t:   size = sizeof(float);   break;
        case DOUBLE_t:  size = sizeof(double);  break;
        case CHAR_t:    size = sizeof(char);    break;
        case INT_t:     size = sizeof(int);     break;
        case BOOL_t:    size = sizeof(bool);    break;
        default: // other types no supported
            size = 0;
            return 0;
        }


        // Multiply by the number of input elements
        if (inputs[0]->ValueType == GRID) {
        	size *= 50625;
        }
        else {
        	size *= inputs[0]->VectorLength;
        }

        // If the inputs are complex values, the total datasize is twice sizeof(VALUE)
        if (inputs[0]->ValueType == VALUE_CMPX) {
            size <<= 1; // times 2 for complex
        }


        // Check if the current buffer is full, and switch to a new one if not
        if ((bufEnd[currUpdateBuf] + size) >= bufSize) {
            // signal thread to start writing previous update Buf
            errCheckMod(osindep_sem_post(&bufFullSem));
            currUpdateBuf++;
            if (currUpdateBuf >= numBufs) currUpdateBuf = 0;
            bufEnd[currUpdateBuf] = 0;
        }

        // Advance the pointer pointing to the beginning of the destination memory for copying
        // Copy the data into the buffer so the Update function can return
        // Note that the buffers are all one contiguous chunk of memory only segmented by the way pointers are handled
        // So, dest = the beginning of the chunk of memory + advance through other buffers ahead of you
        //            + advance past the samples already written into this buffer
        void *dest = (void*)((char*)buf_h + currUpdateBuf*bufSize + bufEnd[currUpdateBuf]);
        if (inputs[0]->MemLoc == CUDA_DEVICE){
            cudaStream_t *cuStream = (cudaStream_t*)cuFlowStream;
            cuCheckM(cudaMemcpyAsync(dest, inputs[0]->Data, size,
                     cudaMemcpyDeviceToHost, *cuStream));
        } else {
            memcpy(dest, inputs[0]->Data, size);
        }

        // Advance the starting location of the buffer for the next iteration
        bufEnd[currUpdateBuf] += size;
        return 0;
    }
    else {
        std::cerr << "[" << ModuleName << "] Update: Thread not running."
                 << std::endl;
        return -1;
    }
}

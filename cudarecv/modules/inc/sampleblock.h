
#ifndef INC__SAMPLEBLOCK_H_
#define INC__SAMPLEBLOCK_H_

#include <cuda_runtime.h>
#include <cstdio>
#include <pthread.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
//#include <semaphore.h>
#include "sem.h"
#include "dsp.h"
#include "module.h"

#define SAMPLE_INPUT_SOURCE_FILE    0
#define SAMPLE_INPUT_SOURCE_SOCKET  1

/** \brief Module that provides a 1ms block of raw samples for processing.
 *
 *  Current implementation only allows 1 SampleBlock object for single
 *  receiver algorithms.
 */
class dsp::SampleBlock : public dsp::Module {

    public:

        SampleBlock(void);

        ~SampleBlock(void);

        /** \brief Updates the output port *Data pointer.
         *
         *  Blocking call till a new block of samples are available. *Data
         *  pointer is updated to point to the new block of samples for
         *  processing.
         *  \return 0 on success. -1 on fail(not running).
         */
        int Update(void* cuFlowStream);

        /** \brief Starts the sampling thread. */
        int Start(void* cuFlowStream);

        /** \brief Stops the sampling thread. */
        int Stop(void);

    private:

        /** \brief Opens input raw data pipe/file and initializes buffers
         *
         *  Expects Filename, SamplingFrequency, and RealTime variables to be
         *  initialized.
         *  \return 0 on success
         */
        int Open(void);

        /** \brief Closes the raw data file. */
        int Close(void);

        /** \brief Number of SampleBlock Instances
         *  Used to keep track of syncing between multiple SampleBlocks. */
        static unsigned char NumInstances;

        /** \brief For multi-thread safe operations with static variables. */
        static pthread_mutex_t mtx;

        /** \brief Used to synchronize between thread and Update. */
        osindep_sem_t samplesAvailSem;
        osindep_sem_t buffersAvailSem;

        static const unsigned int FilenameCapacity = 1024;
        char Filename[FilenameCapacity] = {0};

        /** \brief Number of buffers to create. */
        signed short NumBlocks;

        const signed short NumBlocksDefault = 32; //4;

        /** \brief Array length of 1 buffer corresponding to 1ms. */
        unsigned short BlockLength;

        double SamplingFrequency;
        
        double SampleLength;

        unsigned char InputSourceType;

        // TCP Sockets
        int PortNo = 1111;
        static const unsigned int HostnameCapacity = 1024;
        char Hostname[HostnameCapacity] = {0};
        struct sockaddr_in ServerAddr;

        /** \brief Index to current buffer that the thread is loading. */
        signed short CurrLoadBlock;

        /** \brief Index to current buffer that memory is being copied to cuda device. */
        signed short CurrLoadCudaBlock;

        /** \brief Index to current buffer that dsp::flow is reading from. */
        signed short CurrProcBlock;

        /** \brief If true, real time input stream,
         *         if false, post process from data file. */
        bool RunLive;
        
        /** \brief Number of bytes read prior to the start of sampling
         *         (eg, if starting DPE from a scalar lock). */
        int StartByte = 0;

        /** \brief File descriptor of file being opened. */
        int fd = -1;
        
        /** \brief The offset the pointer is at from the beginning of the file. */
        int skipCheck = -1;

        /** \brief Array of buffers. 1 buffer is an array too. */
        void **BlockArr;

        /** \brief Array of buffers. 1 buffer is an array too. */
        void **BlockArrCuda;

        /** \brief Cuda Stream for Async Memory Transfers. */
        cudaStream_t cuStream;

        /** \brief Thread object that runs GetSamplesThread. */
        pthread_t thread;

        /** \brief Flag to keep the thread loop running. */
        volatile bool KeepRunning = 0;
        volatile bool FirstUpdate = true;

        /** \brief An abstraction to use pthreads with C++ objects. */
        static void * SampleThreadEntry(void * thisObj);

        /** \brief Function that is to be implemented as a thread. */
        void GetSamplesThread(void);


        // Debugging
        int test;

};

#endif

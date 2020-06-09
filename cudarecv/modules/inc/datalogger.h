
#ifndef INC__DATALOGGER_H_
#define INC__DATALOGGER_H_

#include <cuda_runtime.h>
#include <cuComplex.h>
#include <stdint.h>
#include <stdio.h>
#include "sem.h"
#include "dsp.h"
#include "module.h"

/** \brief Correlator module. */
class dsp::DataLogger : public dsp::Module {

    public:

        DataLogger();

        DataLogger(std::string Name);

        ~DataLogger();

        /** \brief Any other initialization before running. */
        int Start(void* cuFlowStream);

        /** \brief Any de-initialization before de-constructor. */
        int Stop(void);

        /** \brief . */
        int Update(void* cuFlowStream);

    protected:

        void ConstructorHelper(void);

        /** \brief Number of channels scalar tracking is currently supporting. */
        unsigned char NumChan;

        cudaStream_t cuStream;

        cudaEvent_t cuEvent;

        unsigned char MaxNumChan_h;

        static const unsigned int FilenameCapacity = 1024;
        char Filename[FilenameCapacity] = "";
        FILE *fd = NULL;

        int Open(void);
        int Close(void);

        size_t bufSize = 1024 * sizeof(float);
        static const unsigned char numBufs = 2;

        // Keep track of how far the buffer filled up, in case it wasn't completely full
        unsigned int bufEnd[numBufs] = {0,0};

        // Pointers to the contiguous chunk of memory of the buffer
        void *buf_d;
        void *buf_h;

        // The index of the buffer currently being written to
        int currUpdateBuf;

        int currThreadBuf;

        bool csv;

        osindep_sem_t bufFullSem;

        /** \brief Thread object that runs GetSamplesThread. */
        pthread_t thread;

        /** \brief Flag to keep the thread loop running. */
        volatile bool KeepRunning = 0;

        /** \brief An abstraction to use pthreads with C++ objects. */
        static void * LoggerThreadEntry(void * thisObj);

        /** \brief Function that is to be implemented as a thread. */
        void LoggerThread(void);

};

#endif

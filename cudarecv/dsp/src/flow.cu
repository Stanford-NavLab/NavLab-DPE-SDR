
#include <cuda_runtime.h>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <fcntl.h>
#include <unistd.h>
#include <sched.h>
#include <pthread.h>
#include <sys/time.h>
#include "dsp.h"
#include "flow.h"
#include "errorhandler.h"

#define FLOW_THREAD_PRIORITY    26
// Note that, in cuCheckHelper, the second parameter "failStatement" is run as failStatement;
// Since this is a macro, the statements "KeepRunning=false; break" are inserted into "failStatement"
// allowing these two lines of code to be executed.
#define cuCheckF(stmt)       cuCheckHelper(stmt, KeepRunning=false; break)

dsp::Flow::~Flow(){} // implementing pure virtual destructor 21109417

int dsp::Flow::LoadFlow(const char *filename) {
    return -1;
}

int dsp::Flow::Start() {
    int i,j;

    std::clog << "[FLOW] Starting." << std::endl;

    cuCheck(cudaStreamCreate(&cuStream));

    for (i=0; i<(int)Mods.size(); i++) {
        j = Mods[i]->Start((void*)&cuStream);
        if (j) {
            std::cerr << "[Flow] Unable to start module "
                      << Mods[i]->GetModuleName() << std::endl;
            for (j=0; j<i; j++)
                Mods[i]->Stop();
            return -1;
        }
    }

    cuCheck(cudaEventCreate(&cuEvent));

    pthread_attr_t attr;
    struct sched_param param;
    if (Mods.size() < 1) {
        std::cerr << "[Flow] Flow not loaded." << std::endl;
        return -1;
    }
    i = pthread_attr_init(&attr);
    if (i) {
        std::cerr << "[Flow] Unable to init thread attr."
                  << std::endl;
        return -1;
    }
    int minPrio = sched_get_priority_min(SCHED_RR);
    int maxPrio = sched_get_priority_max(SCHED_RR);
    /* Posix min requirements is 32 priority levels, Although linux provides
       99 levels, good to check for min/max available levels. */
    param.sched_priority = (maxPrio - minPrio) *
                           FLOW_THREAD_PRIORITY / 32 + minPrio;
    i = pthread_attr_setschedpolicy(&attr, SCHED_RR);
    if (i) {
        std::cerr << "[Flow] Unable to set sample thread real-time."
                  << std::endl;
        return -1;
    }
    i = pthread_attr_setschedparam(&attr, &param);
    if (i) {
        std::cerr << "[Flow] Unable to set sample thread priority."
                  << std::endl;
        return -1;
    }
    i = pthread_create(&thread, &attr, FlowThreadEntry, this);
    if (i) {
        std::cerr << "[Flow] Unable to create thread."
                  << std::endl;
        return -1;
    }

    std::clog << "[FLOW] Started." << std::endl;
    return 0;
}

int dsp::Flow::Stop(void) {
    if (KeepRunning || FlowDone) {
        std::clog << "[Flow] Stopping Flow." << std::endl;
        KeepRunning = false;
        return pthread_join(thread, NULL);
    } else {
        std::clog << "[Flow] Stop: Flow wasn't running." << std::endl;
        return 0;
    }
}

void * dsp::Flow::FlowThreadEntry(void * thisObj) {
    ((dsp::Flow *)thisObj)->FlowThread();
    return NULL;
}

void dsp::Flow::FlowThread() {
    int i, j;
    unsigned long runCount = 0;
    struct timeval start, stop, iterStart, iterStop;
    unsigned long long totalIterDuration = 0; //microseconds
    unsigned int avgIterDuration = 0; //microseconds/
    unsigned char numMax = 40;
    unsigned int maxIterDuration[numMax];
    unsigned long maxCount[numMax];
    unsigned int minCount;
    unsigned int minIterDuration = (unsigned int)-1;
    unsigned int iterDuration;
    unsigned long microsecs;
    double seconds;
    for (i=0; i<numMax; i++) maxIterDuration[i] = 0;
    j = gettimeofday(&start, NULL);
    KeepRunning = true;
    while (KeepRunning) {
        for (i=0; i<(int)Mods.size(); i++) {
            j = Mods[i]->Update((void*)&cuStream);
            if (j) {
                std::cerr << "[Flow] " << Mods[i]->GetModuleName()
                          << "->Update() Failed. \nStopping Flow."
                          << std::endl;
                KeepRunning = false;
                break;
            }
            if (i == 0){
                // start timer after SampleBlock returned
                j = gettimeofday(&iterStart, NULL);
            }
        }
        cuCheckF(cudaStreamSynchronize(cuStream));
        j = gettimeofday(&iterStop, NULL);
        if ((iterStop.tv_sec - iterStart.tv_sec) > 0) {
            iterDuration = 1000000 + iterStop.tv_usec - iterStart.tv_usec
                        + (iterStop.tv_sec - iterStart.tv_sec - 1) * 1000000;
        } else {
            iterDuration = (iterStop.tv_usec - iterStart.tv_usec);
        }
        runCount++;
        for (i=0; i<numMax; i++) {
            if (iterDuration > maxIterDuration[i]) {
                for (j=numMax-1; j>i; j--){
                    maxIterDuration[j] = maxIterDuration[j-1];
                    maxCount[j] = maxCount[j-1];
                }
                maxIterDuration[i] = iterDuration;
                maxCount[i] = runCount;
                break;
            }
        }
        if (iterDuration < minIterDuration) {
            minIterDuration = iterDuration;
            minCount = runCount;
        }
        totalIterDuration += iterDuration;
        //if ((runCount % 2000) == 0)
        //    std::clog << "[Flow] runCount = " << runCount << std::endl;
    }

    std::clog << "[Flow] Flow Stopped." << std::endl;
    std::clog << "[Flow] Stopping Modules." << std::endl;
    for (i=0; i<(int)Mods.size(); i++) {
        Mods[i]->Stop();
    }

    std::clog << "[Flow] runCount = " << runCount << std::endl;
    avgIterDuration = totalIterDuration / runCount;
    std::clog << "[Flow] Average 1ms block duration = " << avgIterDuration
              << " us" << std::endl;
    for (i=0; i<1; i++) //numMax; i++)
        std::clog << "[Flow] Max[" << i << "] 1ms block duration = "
                  << maxIterDuration[i]
                  << " us, run count = " << maxCount[i] << std::endl;
    std::clog << "[Flow] Min 1ms block duration = " << minIterDuration
              << " us, run count = " << minCount << std::endl;
    j = gettimeofday(&stop, NULL);
    if (stop.tv_usec < start.tv_usec) {
        seconds = (double)(stop.tv_sec - start.tv_sec - 1);
        microsecs = 1000000 + stop.tv_usec - start.tv_usec;
    } else {
        seconds = (double)(stop.tv_sec - start.tv_sec);
        microsecs = stop.tv_usec - start.tv_usec;
    }
    seconds += (double)microsecs / 1000000.0;
    std::clog << "[Flow] Total time = " << seconds << " seconds." << std::endl;

    cuCheckV(cudaEventDestroy(cuEvent));
    cuCheckV(cudaStreamDestroy(cuStream));

    FlowDone = true;
}

bool dsp::Flow::CheckFlowState(void) {
	return FlowDone;
}

int dsp::Flow::GetModID(const std::string &ModName)const {
    for (int i=0; i<(int)Mods.size(); i++){
        if (ModName.compare(Mods[i]->GetModuleName()) == 0)
            return i;
    }
    std::cerr << "[FLOW] GetModID: Invalid Module Name: " << ModName << std::endl;
    return -1;
}

int dsp::Flow::SetModParam(const std::string &ModName, const std::string key, const int val) {
    int ModID = GetModID(ModName);
    if (ModID < 0) return -1;
    return Mods[ModID]->SetParam(key, val);
}

int dsp::Flow::SetModParam(const std::string &ModName, const std::string key, const char val) {
    int ModID = GetModID(ModName);
    if (ModID < 0) return -1;
    return Mods[ModID]->SetParam(key, val);
}

int dsp::Flow::SetModParam(const std::string &ModName, const std::string key, const float val) {
    int ModID = GetModID(ModName);
    if (ModID < 0) return -1;
    return Mods[ModID]->SetParam(key, val);
}

int dsp::Flow::SetModParam(const std::string &ModName, const std::string key, const double val) {
    int ModID = GetModID(ModName);
    if (ModID < 0) return -1;
    return Mods[ModID]->SetParam(key, val);
}

int dsp::Flow::SetModParam(const std::string &ModName, const std::string key, const bool val) {
    int ModID = GetModID(ModName);
    if (ModID < 0) return -1;
    return Mods[ModID]->SetParam(key, val);
}

int dsp::Flow::SetModParam(const std::string &ModName, const std::string key, const char *str) {
    int ModID = GetModID(ModName);
    if (ModID < 0) return -1;
    return Mods[ModID]->SetParam(key, str);
}

int dsp::Flow::GetModParam(const std::string &ModName, const std::string key, int *val) {
    int ModID = GetModID(ModName);
    if (ModID < 0) return -1;
    return Mods[ModID]->GetParam(key, val);
}

int dsp::Flow::GetModParam(const std::string &ModName, const std::string key, float *val) {
    int ModID = GetModID(ModName);
    if (ModID < 0) return -1;
    return Mods[ModID]->GetParam(key, val);
}

int dsp::Flow::GetModParam(const std::string &ModName, const std::string key, double *val) {
    int ModID = GetModID(ModName);
    if (ModID < 0) return -1;
    return Mods[ModID]->GetParam(key, val);
}

int dsp::Flow::GetModParam(const std::string &ModName, const std::string key, bool *val) {
    int ModID = GetModID(ModName);
    if (ModID < 0) return -1;
    return Mods[ModID]->GetParam(key, val);
}

int dsp::Flow::GetModParam(const std::string &ModName, const std::string key, char *str,
                           const unsigned int capacity) {
    int ModID = GetModID(ModName);
    if (ModID < 0) return -1;
    return Mods[ModID]->GetParam(key, str, capacity);
}

int dsp::Flow::ConnectPort(const std::string &srcModName, const char *srcPortName,
                           const std::string &dstModName, const char *dstPortName)
{
    int srcModID, srcPortID, dstModID, dstPortID;
    srcModID = GetModID(srcModName);
    if (srcModID < 0) return -1;
    srcPortID = Mods[srcModID]->GetOutputID(srcPortName);
    if (srcPortID < 0) return -1;
    dstModID = GetModID(dstModName);
    if (dstModID < 0) return -1;
    dstPortID = Mods[dstModID]->GetInputID(dstPortName);
    if (dstPortID < 0) return -1;
    return ConnectPort(srcModID, srcPortID, dstModID, dstPortID);
}

int dsp::Flow::ConnectPort(int srcModID, int srcPortID,
                           int dstModID, int dstPortID)
{
    int i;
    dsp::Port *temp;

    if ((srcModID >= (int)Mods.size()) || (dstModID >= (int)Mods.size())) {
        std::cerr << "[Flow] Module IDs do not exist" << std::endl;
        return -1;
    }

    i = Mods[srcModID]->GetOutput(srcPortID, &temp);
    if (i) {
        std::cerr << "[Flow] Unable to get output from ["
                  << Mods[srcModID]->GetModuleName() << "]" << std::endl;
        return -1;
    }

    i = Mods[dstModID]->SetInput(dstPortID, temp);
    if (i) {
        std::cerr << "[Flow] Unable to set input to ["
                  << Mods[dstModID]->GetModuleName() << "]" << std::endl;
        return -1;
    }

    /** \todo Implement some kind of graph data structure to save connections
     *        to allow for generating a graphical block diagram in future.
     */

    return 0;
}

int
dsp::Flow::GetOutput (const std::string &modName,const std::string &portName,dsp::Port **out)const{
    const int modID = GetModID(modName);
    if (modID < 0) return -1;

    char* portname_c = new char[portName.size()+1];
    strcpy(portname_c,portName.c_str());
    const int portID = Mods[modID]->GetOutputID(portname_c);
    delete[] portname_c;

    if (portID < 0) return -1;
    return Mods[modID]->GetOutput(portID,out);
}


#ifndef INC__FLOW_H_
#define INC__FLOW_H_

#include <cuda_runtime.h>
#include <vector>
#include <string>
#include <pthread.h>
#include "dsp.h"
#include "module.h"

namespace dsp {

/** \brief Implements the block/flow diagram of the receiver algorithm. */
class Flow {
    friend class dsp::FlowMgr;

    public:

        /** \brief Loads the flow diagram from file.
         *  \param filename Path to the text file to load from.
         *  \return 0 on success. -1 on fail.
         */
        virtual int LoadFlow(const char *filename);

        /** \brief Starts the Flow thread. */
        virtual int Start(void);

        /** \brief Stops the Flow thread. */
        virtual int Stop(void);

        /** \brief Queries the state of the Flow thread. */
        virtual bool CheckFlowState(void);

        /** \brief Virtual destructor */
        virtual ~Flow() = 0; // "-Wdelete-non-virtual-dtor"

        int     GetOutput (const std::string&,const std::string&,dsp::Port**)const;

        // Abuse of power so that main can run multiple flows -- these should really be protected
        //int SetModParam(const std::string &ModName, const std::string key, Param *p);
        int SetModParam(const std::string &ModName, const std::string key, const int val);
        int SetModParam(const std::string &ModName, const std::string key, const char val);
        int SetModParam(const std::string &ModName, const std::string key, const float val);
        int SetModParam(const std::string &ModName, const std::string key, const double val);
        int SetModParam(const std::string &ModName, const std::string key, const bool val);
        int SetModParam(const std::string &ModName, const std::string key, const char *str);



    protected:

        /** \brief To be used as an array of Modules used in this flow. */
        std::vector<dsp::Module*> Mods;

        pthread_t thread;

        /** \brief Flag to keep the thread loop running. */
        volatile bool KeepRunning = 0;

        /** \brief Flag to signify when the flow has ended. */
        volatile bool FlowDone = false;

        /** \brief An abstraction to use pthreads with C++ objects. */
        static void * FlowThreadEntry(void * thisObj);

        /** \brief The actual thread that will run. */
        void FlowThread(void);

        cudaStream_t cuStream;
        cudaEvent_t cuEvent;

        int GetModID(const std::string &ModName)const;

        //int GetModParam(const std::string &ModName, const std::string key, Param *p);
        int GetModParam(const std::string &ModName, const std::string key, int *val);
        int GetModParam(const std::string &ModName, const std::string key, float *val);
        int GetModParam(const std::string &ModName, const std::string key, double *val);
        int GetModParam(const std::string &ModName, const std::string key, bool *val);
        int GetModParam(const std::string &ModName, const std::string key, char *str, const unsigned int capacity);

        /** \brief Connect output port of module to input port of module
         *
         *  \param srcModName Module to get output from
         *  \param srcPortName Output port
         *  \param dstModName Module to connect port into
         *  \param dstPortName Input port
         *  \return 0 on success, -1 otherwise
         */
        int ConnectPort(const std::string &srcModName, const char *srcPortName,
                        const std::string &dstModName, const char *dstPortName);

        /** \brief Connect output port of module to input port of module
         *
         *  \param srcModID Module number of which to get output from
         *  \param srcPortID Output port number
         *  \param dstModID Module number to connect port into
         *  \param dstPortID Input port number
         *  \return 0 on success, -1 otherwise
         */
        int ConnectPort(int srcModID, int srcPortID,
                        int dstModID, int dstPortID);

};

}

#endif

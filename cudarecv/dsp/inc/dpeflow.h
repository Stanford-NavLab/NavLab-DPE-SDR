
#ifndef INC_DPEFLOW_H_
#define INC_DPEFLOW_H_

#include <pthread.h>
#include "dsp.h"
#include "flow.h"

class dsp::DPEFlow : public dsp::Flow {

    public:
        virtual ~DPEFlow(){}

        /** \brief Loads hardcoded scalar tracking.
         *  \param filename dummy input arg.
         *  \return 0 on success. -1 on fail.
         */
        int LoadFlow(const char *filename);

};

#endif


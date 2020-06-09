
#ifndef INC__MODULE_H_
#define INC__MODULE_H_

#include <string>
#include <map>
#include <iterator>
#include "dsp.h"

namespace dsp {

/** \brief Base class that all modules should inherit from. */
class Module {

    public:

        /** \brief Processes information for that time instance.
         *
         *  Processes input(s) and produces output(s) for that time instance.
         *  For example a 1ms block.
         *  \return 0 if no errors.
         */
        virtual int Update(void* cuFlowStream) = 0;

        int GetInputID(const char *inName);

        int GetOutputID(const char *outName);

        int SetInput(const dsp::Port &in);

        /** \brief "Connects" the input port.
         *
         *  Assigns an existing port as an input to this module.
         *  \param id Input port number to connect to.
         *  \param *in Pointer to the port.
         *  \return 0 on success, -1 otherwise.
         */
        int SetInput(unsigned char id, dsp::Port *in);

        /** \brief "Connects" the output port
         *
         *  Obtains the pointer to the output port of the module.
         *  \param id Output port number.
         *  \param **out Pointer to store the pointer to.
         *  \return 0 on success, -1 otherwise.
         */
        int GetOutput(unsigned char id, dsp::Port **out);

        /** \brief Get name of the module/ child class. */
        std::string GetModuleName(void);

        int SetParam(const std::string key, Param *p);
        int SetParam(const std::string key, const int val);
        int SetParam(const std::string key, const char val);
        int SetParam(const std::string key, const float val);
        int SetParam(const std::string key, const double val);
        int SetParam(const std::string key, const bool val);
        int SetParam(const std::string key, const char *str);

        int GetParam(const std::string key, Param *p);
        int GetParam(const std::string key, int *val);
        int GetParam(const std::string key, float *val);
        int GetParam(const std::string key, double *val);
        int GetParam(const std::string key, bool *val);
        int GetParam(const std::string key, char *str, const unsigned int capacity);

        /** \brief Initializes/Starts the module.
         *
         *  Does nothing unless child class implements it.
         *  Meant as a function to allow any initialization that required for
         *  any particular module. It could also be used to spawn a thread if
         *  needed by that module.
         *  \return 0 on success.
         */
        virtual int Start(void* cuFlowStream);

        /** \brief Stops the module.
         *
         *  Does nothing unless implemented by the child module. Child modules
         *  that implement Start should have implemented a corresponding Stop
         *  too. This could include wrapping up any internal states or stopping
         *  a module specific thread.
         *  \return 0 on success.
         */
        virtual int Stop(void);

    protected:

        /** \brief Name of module/ child class. */
        std::string ModuleName;

        std::map<std::string, Param> Params;

        int InsertParam(const std::string key, void *ptr, DataType_t dtype,
                        unsigned int capacity, unsigned int size);

        std::map<std::string, dsp::Port> Inputs;

        std::map<std::string, dsp::Port> Outputs;

        int InsertOutput(const std::string key, DataType_t dtype, ValueType_t vType,
                         MemLoc_t MemLoc, unsigned short VectorLength, void *Data,
                         int AuxValue);

        /** \brief Number of input ports this module is accepting. */
        unsigned char NumInputs = 0;

        /** \brief Number of output ports this module is producing. */
        unsigned char NumOutputs = 0;

        int AllocateInputs(unsigned char num);

        int AllocateOutputs(unsigned char num);

        int ConfigExpectedInput(const unsigned char id, const char *name,
                                const DataType_t dtype, const ValueType_t vType,
                                const unsigned short VectorLength);

        int ConfigOutput(const unsigned char id, const char *name,
                         const DataType_t dtype, const ValueType_t vType,
                         const MemLoc_t MemLoc, const unsigned short VectorLength,
                         void *Data, const int AuxValue);

        int UpdateOutput(const unsigned char id, const unsigned short VectorLength,
                         void *Data, const int AuxValue);

        /** \brief An Array of expected inputs.
         *
         *  This array provides information on what kind of inputs the module
         *  is accepting for each corresponding port number. This allows
         *  SetInput() to provide error checking for appropriate port
         *  connections.
         */
        dsp::ExpectedPort *expectedInputs = NULL;

        /** \brief Array of input port pointers. */
        dsp::Port **inputs = NULL;

        /** \brief Array of outputs. */
        dsp::Port *outputs = NULL;

        virtual ~Module () = 0;

};

}

#endif

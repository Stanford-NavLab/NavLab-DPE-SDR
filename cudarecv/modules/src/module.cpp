
#include <iostream>
#include <cstring>
#include "module.h"

int dsp::Module::InsertOutput(const std::string key, DataType_t dtype, ValueType_t vType,
                         MemLoc_t MemLoc, unsigned short VectorLength, void *Data,
                         int AuxValue) {
    /*std::pair<std::map<std::string, dsp::Port>::iterator,bool> ret;
    dsp::Port p;
    ret = Outputs.insert(std::pair<std::string,dsp::Port>(key,p));
    if (ret.second==false) {
        std::cerr << "[" << ModuleName
                  << "] InsertPort: Parameter with key \""
                  << key << "\" already exists." << std::endl;
        return -1;
    }*/
    return -1;
}

int dsp::Module::InsertParam(const std::string key, void *ptr, DataType_t dtype,
                             unsigned int capacity, unsigned int size) {
    std::pair<std::map<std::string,Param>::iterator,bool> ret;
    Param p;
    p.Ptr = ptr;
    p.Datatype = dtype;
    p.Capacity = capacity;
    p.Size = size;
    ret = Params.insert(std::pair<std::string,Param>(key,p));
    if (ret.second==false) {
        std::cerr << "[" << ModuleName
                  << "] InsertParam: Parameter with key \""
                  << key << "\" already exists." << std::endl;
        return -1;
    }
    return 0;
}


int dsp::Module::SetParam(const std::string key, Param *p) {
    std::map<std::string,Param>::iterator it = Params.find(key);
    if (it == Params.end()) {
        std::cerr << "[" << ModuleName
                  << "] SetParam: Parameter with key \""
                  << key << "\" does not exist." << std::endl;
        return -1;
    }
    if (it->second.Datatype != p->Datatype) {
        std::cerr << "[" << ModuleName
                  << "] SetParam: Parameter with key \""
                  << key << "\" datatype does not match." << std::endl;
        return -1;
    }
    if (it->second.Capacity < p->Size) {
        std::cerr << "[" << ModuleName
                  << "] SetParam: Parameter with key \""
                  << key << "\" has size larger than capacity." << std::endl;
        return -1;
    }
    it->second.Size = p->Size;
    memcpy(it->second.Ptr, p->Ptr, p->Size);
    return 0;
}

// The following member functions catch the type of the variable passed in
// by the user command and convert it into a type param for the above
// member function to process
int dsp::Module::SetParam(const std::string key, const int val) {
    Param p;
    p.Ptr = (void*)&val;
    p.Datatype = INT_t;
    p.Size = sizeof(int);
    return SetParam(key, &p);
}

int dsp::Module::SetParam(const std::string key, const char val) {
    Param p;
    p.Ptr = (void*)&val;
    p.Datatype = CHAR_t;
    p.Size = sizeof(char);
    return SetParam(key, &p);
}

int dsp::Module::SetParam(const std::string key, const float val) {
    Param p;
    p.Ptr = (void*)&val;
    p.Datatype = FLOAT_t;
    p.Size = sizeof(float);
    return SetParam(key, &p);
}

int dsp::Module::SetParam(const std::string key, const double val) {
    Param p;
    p.Ptr = (void*)&val;
    p.Datatype = DOUBLE_t;
    p.Size = sizeof(double);
    return SetParam(key, &p);
}

int dsp::Module::SetParam(const std::string key, const bool val) {
    Param p;
    p.Ptr = (void*)&val;
    p.Datatype = BOOL_t;
    p.Size = sizeof(bool);
    return SetParam(key, &p);
}

int dsp::Module::SetParam(const std::string key, const char *str) {
    Param p;
    p.Ptr = (void*)str;
    p.Datatype = CHAR_t;
    p.Size = sizeof(char) * (strlen(str) + 1);  // copy null character too
    return SetParam(key, &p);
}

int dsp::Module::GetParam(const std::string key, Param *p) {
    std::map<std::string,Param>::iterator it = Params.find(key);
    if (it == Params.end()) {
        std::cerr << "[" << ModuleName
                  << "] GetParam: Parameter with key \""
                  << key << "\" does not exist." << std::endl;
        return -1;
    }
    if (p->Capacity < it->second.Size) {
        std::cerr << "[" << ModuleName
                  << "] GetParam: Parameter with key \""
                  << key << "\" has size larger than capacity." << std::endl;
        return -1;
    }
    p->Datatype = it->second.Datatype;
    p->Size = it->second.Size;
    memcpy(p->Ptr, it->second.Ptr, p->Size);
    return 0;
}

int dsp::Module::GetParam(const std::string key, int *val) {
    Param p;
    p.Ptr = (void*)val;
    p.Datatype = INT_t;
    p.Capacity = sizeof(int);
    return GetParam(key, &p);
}

int dsp::Module::GetParam(const std::string key, float *val) {
    Param p;
    p.Ptr = (void*)val;
    p.Datatype = FLOAT_t;
    p.Capacity = sizeof(float);
    return GetParam(key, &p);
}

int dsp::Module::GetParam(const std::string key, double *val) {
    Param p;
    p.Ptr = (void*)val;
    p.Datatype = DOUBLE_t;
    p.Capacity = sizeof(double);
    return GetParam(key, &p);
}

int dsp::Module::GetParam(const std::string key, bool *val) {
    Param p;
    p.Ptr = (void*)val;
    p.Datatype = BOOL_t;
    p.Capacity = sizeof(bool);
    return GetParam(key, &p);
}

int dsp::Module::GetParam(const std::string key, char *str,
                          const unsigned int capacity) {
    int ret;
    Param p;
    p.Ptr = (void*)str;
    p.Datatype = CHAR_t;
    p.Capacity = sizeof(char) * capacity;
    ret =  GetParam(key, &p);
    str[capacity-1] = 0; //null terminate if overflow
    return ret;
}

int dsp::Module::AllocateInputs(unsigned char num) {
    if ((NumInputs != 0) || (expectedInputs != NULL) || (inputs != NULL)) {
        std::cerr << "[" << ModuleName
                  << "] AllocateInputs: Inputs already allocated" << std::endl;
        return -1;
    }
    NumInputs = num;
    if (num > 0) {
        expectedInputs = new dsp::ExpectedPort[num];
        inputs = new dsp::Port *[num];
    }
    return 0;
}

int dsp::Module::
ConfigExpectedInput(const unsigned char id, const char *name, const DataType_t dtype,
                    const ValueType_t vType, const unsigned short VectorLength) {
    if (id >= NumInputs) {
        std::cerr << "[" << ModuleName
                  << "] ConfigExpectedInput: id " << name << " out of range." << std::endl;
        return -1;
    }
    std::strcpy(expectedInputs[id].Name, name);
    expectedInputs[id].ValueType = vType;
    expectedInputs[id].VectorLength = VectorLength;
    expectedInputs[id].Datatype = dtype;
    return 0;
}

int dsp::Module::
ConfigOutput(const unsigned char id, const char *name, const DataType_t dtype,
             const ValueType_t vType, const MemLoc_t MemLoc,
             const unsigned short VectorLength, void *Data, const int AuxValue) {
    if (id >= NumOutputs) {
        std::cerr << "[" << ModuleName
                  << "] ConfigOutput: id out of range." << std::endl;
        return -1;
    }
    std::strcpy(outputs[id].Name, name);
    outputs[id].Datatype = dtype;
    outputs[id].Exponent = 0;
    outputs[id].ValueType = vType;
    outputs[id].VectorLength = VectorLength;
    outputs[id].AuxValue = AuxValue;
    outputs[id].Data = Data;
    outputs[id].MemLoc = MemLoc;
    return 0;
}

int dsp::Module::
UpdateOutput(const unsigned char id, const unsigned short VectorLength,
             void *Data, const int AuxValue) {
    if (id >= NumOutputs) {
        std::cerr << "[" << ModuleName
                  << "] UpdateOutput: id out of range." << std::endl;
        return -1;
    }
    outputs[id].VectorLength = VectorLength;
    outputs[id].AuxValue = AuxValue;
    outputs[id].Data = Data;
    return 0;
}

int dsp::Module::AllocateOutputs(unsigned char num) {
    if ((NumOutputs != 0) || (outputs != NULL)) {
        std::cerr << "[" << ModuleName
                  << "] AllocateOutputs: Outputs already allocated" << std::endl;
        return -1;
    }
    NumOutputs = num;
    if (num > 0) {
        outputs = new dsp::Port[num];
    }
    return 0;
}

int dsp::Module::GetInputID(const char *inName) {
    for (unsigned char i=0; i<NumInputs; i++){
        if (std::strcmp(inName, expectedInputs[i].Name) == 0)
            return (int)i;
    }
    std::cerr << "[" << ModuleName << "] GetInputID: " << inName
              << " does not exist." << std::endl;
    return -1;
}

int dsp::Module::GetOutputID(const char *outName) {
    for (unsigned char i=0; i<NumOutputs; i++){
        if (std::strcmp(outName, outputs[i].Name) == 0)
            return (int)i;
    }
    std::cerr << "[" << ModuleName << "] GetOutputID: " << outName
              << " does not exist." << std::endl;
    return -1;
}

dsp::Module::~Module(){} // stackoverflow 2555033

int dsp::Module::SetInput(const dsp::Port &in) {
    /** \todo Implement intelligent input port matching. */
    std::cerr << "[Module] Function not implemented." << std::endl;
    return -1;
}

int dsp::Module::SetInput(unsigned char id, dsp::Port *in) {
    if (id < NumInputs) {
        if ((expectedInputs[id].ValueType == in->ValueType) ||
            (expectedInputs[id].ValueType == ANY)) {
            if (((expectedInputs[id].VectorLength == in->VectorLength) &&
                (expectedInputs[id].VectorLength > 0)) ||
                (expectedInputs[id].VectorLength == VECTORLENGTH_ANY)) {
                // Valid Input
                inputs[id] = in;
                return 0;
            } else {
                std::cerr << "[" << ModuleName <<
                            "] Invalid input vector length." << std::endl <<
                            "Input size: " << in->VectorLength << std::endl <<
                            "Expected size: " << expectedInputs[id].VectorLength << std::endl;
                return -1;
            }
        } else {
            std::cerr << "[" << ModuleName <<
                        "] Invalid input value type." << std::endl;
            return -1;
        }
    } else {
        std::cerr << "[" << ModuleName << "] No such input." << std::endl;
        return -1;
    }
}

int dsp::Module::GetOutput(unsigned char id, dsp::Port **out) {
    if (id < NumOutputs) {
        *out = &outputs[id];
        return 0;
    } else {
        std::cerr << "[" << ModuleName << "] No such output." << std::endl;
        return -1;
    }
}

std::string dsp::Module::GetModuleName() {
    return ModuleName;
}


int dsp::Module::Start(void* cuFlowStream) {
    return 0;
}

int dsp::Module::Stop() {
    return 0;
}

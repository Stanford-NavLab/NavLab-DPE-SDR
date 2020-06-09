
#ifndef _FLOWMGR_H
#define _FLOWMGR_H

#include "auxil.h"
#include "dsp.h"
#include "flow.h"


#include "dpeflow.h"

#include <vector>
#include <string>
#include <map>
#include <utility>
#include <signal.h>

#define NPOS size_t(-1)

#define FlowFnt(T,CONST)            \
int T(const std::string& key)CONST{ \
    size_t idx = getFlowIdx (key);  \
    return (idx == NPOS)?-1:T(idx); \
}

#define RegFlow(TYPE,CNAME,M,NAME,DES)          \
class CNAME: public Registrar{                  \
public:                                         \
    CNAME(){}                                   \
    ~CNAME(){}                                  \
    std::string get_name ()const{return NAME;}  \
    std::string get_desc ()const{return DES;}   \
    bool match(const std::string& cmd)const{    \
        return (auxil::strcmpi(NAME,cmd,M) == 0); \
    }                                           \
    Flow* createFlow()const{                    \
        return new dsp::TYPE;                   \
    }                                           \
}

class dsp::FlowMgr {
    using FlowIdent = std::pair<dsp::Flow*,std::string>;
    std::vector<FlowIdent> _flowlist;
    std::map<std::string,size_t> _alias;

    FlowMgr             (const FlowMgr&);
    Flow*   getFlowPtr  (const size_t& idx)const;
    int     loadFlow    (const size_t& , const char* = nullptr) const;

    int     setParam    (const size_t&,const std::string&,const std::string&,const std::string&)const;
    int     listOutput  (const size_t&,const std::string&,const std::string&)const;

    int     destroyFlow (const std::size_t&);
    int     startFlow   (const std::size_t&)const;
    int     stopFlow    (const std::size_t&)const;

public:
    FlowMgr();
    ~FlowMgr();
    int createFlow      (const std::string&, const std::string& = "");

    int     loadFlow    (const std::string& key, const char* = nullptr) const;
    size_t  getFlowIdx  (const std::string&) const;
    void    addAlias    (const std::string& key,const size_t& value);
    void    listAlias   () const;
    void    listFlow    () const;
    void    flowType    () const;
    int     setParam    (const std::string&,const std::string&,const std::string&,const std::string&)const;
    int     listOutput  (const std::string&,const std::string&,const std::string&)const;

    FlowFnt (destroyFlow,     );
    FlowFnt (startFlow,  const);
    FlowFnt (stopFlow,   const);
    size_t  EmergencyStop()const;


private:
    class Registrar;
    std::vector<Registrar*> _regis;

    class Registrar {
    protected:
        Registrar(){}

    public:
        virtual ~Registrar()                              = 0;
        virtual Flow* createFlow()                  const = 0;
        virtual std::string get_name()              const = 0;
        virtual std::string get_desc()              const = 0;
        virtual bool match(const std::string& cmd)  const = 0;
    };

    RegFlow(DPEFlow,    DPE         , 3 ,"DPE"          ,"Direct Position Estimation Flow");
};





#endif

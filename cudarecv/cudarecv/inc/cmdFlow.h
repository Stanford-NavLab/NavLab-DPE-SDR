
#ifndef CMD_FLOW_H
#define CMD_FLOW_H

#include "console.h"
#include "cmdExec.h"
#include "flowmgr.h"
#include <string>
#include <iostream>
#include <iomanip>

#define CmdFlowClass(T,N,US,HE)                             \
namespace console {                                         \
    class T;                                                \
};                                                          \
class console::T : public console::CmdExec{                 \
    dsp::FlowMgr*const _fmgr;                               \
public:                                                     \
    T(console::CmdParser*const cm,dsp::FlowMgr*const fm):   \
        CmdExec(cm),_fmgr(fm){}                             \
    ~T() {}                                                 \
    console::CmdExecStatus exec(const std::string& option); \
    void usage(std::ostream& os) const{                     \
        os << "Usage: " << N << ' ' << US << std::endl;     \
    }                                                       \
    void help() const{                                      \
        std::cout << std::setw(15) << std::left << N        \
              << HE << std::endl;                           \
    }                                                       \
}

namespace console{
    bool initFlowCmd(CmdParser*const,dsp::FlowMgr*const);
};

CmdFlowClass(NewFlow,"NEWFlow","<(flow type)> [(flow alias)]","Create a new flow and optionally assign alias");
CmdFlowClass(DelFlow,"DELFlow","[(flow ident #1) ...]","Delete flow(s) using alias/index");
CmdFlowClass(StartFlow,"STARTFlow","[(flow ident #1) ...]","Start flow(s) using alias/index");
CmdFlowClass(StopFlow,"STOPFlow","[(flow ident #1) ...]","Stop flow(s) using alias/index");
CmdFlowClass(LoadFlow,"LOADFlow","<(flow ident)> [(file path)]","Load a flow using alias/index");
CmdFlowClass(AddAlias,"ADDAlias","<(flow ident)> [(new alias #1) ...]","Add new alias to flow of alias/index");
CmdFlowClass(ActiveAlias,"ACTAlias","","List all aliases-index pairs");
CmdFlowClass(ActiveFlow ,"ACTFlow","","List all active flows");
CmdFlowClass(SetParam ,"SETParam","<(flow ident)> <(module name)> <(parameter name)> <(parameter value)>","Set module parameter");
// TODO: show status; selectively shown according to type.
CmdFlowClass(PrintPort, "PRINTport","<(flow ident)> <(module name)> <(port name)>","Show the output from a specific port");
CmdFlowClass(FlowType, "LSFlow" ,"","Show all flow types");

#endif

/****************************************************************************
  FileName     [ cmdCommon.h ]
  PackageName  [ cmd ]
  Synopsis     [ Define classes for common commands ]
  Author       [ Chung-Yang (Ric) Huang ]
  Copyright    [ Copyleft(c) 2007-2013 LaDs(III), GIEE, NTU, Taiwan ]
****************************************************************************/
#ifndef CMD_COMMON_H
#define CMD_COMMON_H

#include "console.h"
#include "cmdExec.h"
#include <string>

#define CmdClass(T)                              \
class console::T : public console::CmdExec{      \
public:                                          \
    T(console::CmdParser*const cm):CmdExec(cm){} \
    ~T() {}                                      \
    console::CmdExecStatus exec(const std::string& option);  \
    void usage(std::ostream& os) const;          \
    void help() const;                           \
}

CmdClass(HelpCmd);
CmdClass(QuitCmd);
CmdClass(HistoryCmd);
CmdClass(DofileCmd);

#endif // CMD_COMMON_H

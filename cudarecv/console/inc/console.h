
#ifndef _CONSOLE_H
#define _CONSOLE_H
#include <iostream>

namespace dsp {
    class FlowMgr;
}

namespace console {
    class CmdExec;
    class CmdParser;
    char getKey(std::istream& is = std::cin);

    enum CmdExecStatus:size_t {
        CMD_EXEC_QUIT  = 0,
        CMD_EXEC_DONE  = 1,
        CMD_EXEC_ERROR = 2,
        CMD_EXEC_NOP   = 3,

        // dummy
        CMD_EXEC_TOT
    };

    enum CmdOptionError:size_t {
        CMD_OPT_MISSING    = 0,
        CMD_OPT_EXTRA      = 1,
        CMD_OPT_ILLEGAL    = 2,
        CMD_OPT_FOPEN_FAIL = 3,

        // dummy
        CMD_OPT_ERROR_TOT
    };

    /* Common commands */
    bool initCommonCmd(CmdParser*const cm);
    class HelpCmd;
    class QuitCmd;
    class HistoryCmd;
    class DofileCmd;

    bool initFlowCmd(CmdParser*const cm,dsp::FlowMgr*const fm);

    class os { // 772355
        std::ostream &_out;
    public:
        os (std::ostream& out):_out(out){}
        template<typename T>
        console::os& operator << (const T& v);
    };
};

#endif

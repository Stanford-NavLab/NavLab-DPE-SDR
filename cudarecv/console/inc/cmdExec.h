
#ifndef _CMD_EXEC_H
#define _CMD_EXEC_H

#include "console.h"
#include <iostream>
#include <string>
#include <vector>

class console::CmdExec {
public:
    CmdExec(console::CmdParser*const cm);
    virtual ~CmdExec();

    virtual CmdExecStatus exec        (const std::string& str) = 0;
    virtual void          usage       (std::ostream& os) const = 0;
    virtual void          help        ()                 const = 0;

    void                  setOptCmd   (const std::string& str);
    bool                  checkOptCmd (const std::string& check) const;
    const std::string&    getOptCmd   () const;

protected:
    bool            lexSingleOption (const std::string&, std::string&,              bool optional = true) const;
    bool            lexOptions      (const std::string&, std::vector<std::string>&, size_t nOpts  = 0)    const;
    CmdExecStatus   errorOption     (CmdOptionError err, const std::string& opt) const;

    CmdParser*const _cmdMgr;


private:
    std::string                _optCmd;
};

#endif


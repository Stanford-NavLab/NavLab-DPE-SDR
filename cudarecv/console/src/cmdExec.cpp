
#include "console.h"
#include "cmdExec.h"
#include "auxil.h"
//using namespace std;
using namespace console;

console::CmdExec::CmdExec(CmdParser*const cm):_cmdMgr(cm)  {}
console::CmdExec::~CmdExec() {}

void
CmdExec::setOptCmd(const std::string& str) {
    _optCmd = str;
}

const std::string&
CmdExec::getOptCmd() const {
    return _optCmd;
}

//----------------------------------------------------------------------
//     Member Function for class CmdExec
//----------------------------------------------------------------------
// Return false if error options found
// "optional" = true if the option is optional XD
// "optional": default = true
//
bool
CmdExec::lexSingleOption
(const std::string& option, std::string& token, bool optional) const {
    const size_t n = auxil::extractTok(option,token);
    if (!optional) {
        if (token.size() == 0) {
            errorOption(CMD_OPT_MISSING, "");
            return false;
        }
    }
    if (n != std::string::npos) {
        errorOption(CMD_OPT_EXTRA, option.substr(n));
        return false;
    } return true;
}

// if nOpts is specified (!= 0), the number of tokens must be exactly = nOpts
// Otherwise, return false.
//
bool
CmdExec::lexOptions
(const std::string& option, std::vector<std::string>& tokens, size_t nOpts) const {

    std::string tok;
    for (size_t n = 0; n != std::string::npos ; tokens.push_back(tok)){
        n = auxil::extractTok(option, tok, n);
        if (tok.empty()) break;
    }

    if (nOpts == 0 || tokens.size() == nOpts) return true;

    if (tokens.size() < nOpts) errorOption(CMD_OPT_MISSING, "");
    else errorOption(CMD_OPT_EXTRA, tokens[nOpts]);
    return false;
}

CmdExecStatus
CmdExec::errorOption(CmdOptionError err, const std::string& opt) const {
    switch (err) {
        case CMD_OPT_MISSING:
            std::cerr << "Missing option";
            if (opt.size()) std::cerr << " after (" << opt << ")";
            std::cerr << "!!" << std::endl;
        break;
        case CMD_OPT_EXTRA:
            std::cerr << "Extra option!! (" << opt << ")" << std::endl;
        break;
        case CMD_OPT_ILLEGAL:
            std::cerr << "Illegal option!! (" << opt << ")" << std::endl;
        break;
        case CMD_OPT_FOPEN_FAIL:
            std::cerr << "Error: cannot open file \"" << opt << "\"!!" << std::endl;
        break;
        default:
            std::cerr << "Unknown option error type!! (" << err << ")" << std::endl;
        exit(-1);
    } return CMD_EXEC_ERROR;
}

// Called by "getCmd()"
// Check if "check" is a matched substring of "_optCmd"...
// if not, return false.
//
// Perform case-insensitive checks
//
bool
CmdExec::checkOptCmd(const std::string& check) const {
    if (check.empty()) return true;
    if (check.size() > _optCmd.size()) return false;
    return (auxil::strcmpi(_optCmd,check,1)==0);
}


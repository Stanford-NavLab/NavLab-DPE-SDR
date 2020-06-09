/****************************************************************************
  FileName     [ cmdCommon.cpp ]
  PackageName  [ cmd ]
  Synopsis     [ Define common commands ]
  Author       [ Chung-Yang (Ric) Huang ]
  Copyright    [ Copyleft(c) 2007-2013 LaDs(III), GIEE, NTU, Taiwan ]
****************************************************************************/
#include <iomanip>
#include <string>
#include <cstdlib>
#include "auxil.h"
#include "cmdCommon.h"
#include "console.h"
#include "cmdExec.h"
#include "cmdParser.h"

using namespace std;
using namespace console;


bool
console::initCommonCmd(CmdParser*const cm) {
    if (!(cm->regCmd("Quit"     , 1, new QuitCmd    (cm)) &&
          cm->regCmd("HIStory"  , 3, new HistoryCmd (cm)) &&
          cm->regCmd("HELp"     , 3, new HelpCmd    (cm)) &&
          cm->regCmd("DOfile"   , 2, new DofileCmd  (cm))
        )) {
        cerr << "Registering \"init\" commands fails... exiting" << endl;
        return false;
    } return true;
}

//----------------------------------------------------------------------
//     HELp [(string cmd)]
//----------------------------------------------------------------------
CmdExecStatus
HelpCmd::exec(const string& option) {
    // check option
    string token;
    if (!CmdExec::lexSingleOption(option, token))
        return CMD_EXEC_ERROR;
    if (token.size()) {
        CmdExec* e = _cmdMgr->getCmd(token);
        if (!e) return CmdExec::errorOption(CMD_OPT_ILLEGAL, token);
        e->usage(cout);
    }
    else
        _cmdMgr->printHelps();
    return CMD_EXEC_DONE;
}

void
HelpCmd::usage(ostream& os) const {
    os << "Usage: HELp [(string cmd)]" << endl;
}

void
HelpCmd::help() const {
    cout << setw(15) << left << "HELp: "
          << "print this help message" << endl;
}

//----------------------------------------------------------------------
//     Quit [-Force]
//----------------------------------------------------------------------
CmdExecStatus
QuitCmd::exec(const string& option) {
    // check option
    string token;
    if (!CmdExec::lexSingleOption(option, token))
        return CMD_EXEC_ERROR;
    if (token.size()&& token.size() <=7) {
        if (auxil::strcmpi("-Forced",token,2) != 0)
            return CmdExec::errorOption(CMD_OPT_ILLEGAL, token);
        else
            return CMD_EXEC_QUIT;  // ready to quit
    }

    cout << "Are you sure to quit (Yes/No)? [Yes] ";
    string ss;
    getline(cin,ss);
    size_t s = ss.find_first_not_of(' ', 0);
    if (s != string::npos) {
        ss = ss.substr(s);
        if (auxil::strcmpi("No", ss, 1) == 0)
            return CMD_EXEC_DONE;  // ready to quit
    }
    return CMD_EXEC_QUIT;      // not yet to quit
}

void
QuitCmd::usage(ostream& os) const {
    os << "Usage: Quit [-Force]" << endl;
}

void
QuitCmd::help() const {
    cout << setw(15) << left << "Quit: "
          << "quit the execution" << endl;
}

//----------------------------------------------------------------------
//     HIStory [(int nPrint)]
//----------------------------------------------------------------------
CmdExecStatus
HistoryCmd::exec(const string& option) {
    // check option
    string token;
    if (!CmdExec::lexSingleOption(option, token))
        return CMD_EXEC_ERROR;
    int nPrint = -1;
    if (!token.empty()) {
        if (token.find_first_not_of("0123456789") != string::npos)
            return CmdExec::errorOption(CMD_OPT_ILLEGAL, token);
        nPrint = strtol(token.c_str(),nullptr,0);
    } _cmdMgr->printHistory(nPrint);

    return CMD_EXEC_DONE;
}

void
HistoryCmd::usage(ostream& os) const {
    os << "Usage: HIStory [(int nPrint)]" << endl;
}

void
HistoryCmd::help() const {
    cout << setw(15) << left << "HIStory: "
          << "print command history" << endl;
}


//----------------------------------------------------------------------
//     DOfile <(string file)>
//----------------------------------------------------------------------
// TODO: You DON'T need to modify this function!
//         But you should modify CmdParser::openDofile(), etc, in order
//         to support the following features.
//
// Supported features
// (1) mcalc> dofile do1
//      mcalc> ...          <== some other commands
//      mcalc> dofile do2 <== there is a "dofile do1" in do2
//      mcalc>
// (2) mcalc> dofile t
//      Error: cannot open file "t"!!
//      mcalc> dofile do <== can open a dofile "do" after failing to open "t"
//      mcalc>
// (3) If a dofile xx contains a line "dofile xx" calling itself,
//      where xx may or may not exist...  (recursive dofiles)
//      (Let the max recursion depth = 1024)
//
CmdExecStatus
DofileCmd::exec(const string& option) {
    // check option
    string token;
    if (!CmdExec::lexSingleOption(option, token, false))
        return CMD_EXEC_ERROR;
    if (!_cmdMgr->openDofile(token))
        return CmdExec::errorOption(CMD_OPT_FOPEN_FAIL, token);
    return CMD_EXEC_DONE;
}

void
DofileCmd::usage(ostream& os) const {
    os << "Usage: DOfile <(string file)>" << endl;
}

void
DofileCmd::help() const {
    cout << setw(15) << left << "DOfile: "
          << "execute the commands in the dofile" << endl;
}

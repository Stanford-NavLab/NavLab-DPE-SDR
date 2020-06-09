/****************************************************************************
  FileName     [ cmdParser.cpp ]
  PackageName  [ cmd ]
  Synopsis     [ Define command parsing member functions for class CmdParser ]
  Author       [ Chung-Yang (Ric) Huang ]
  Copyright    [ Copyleft(c) 2007-2013 LaDs(III), GIEE, NTU, Taiwan ]
****************************************************************************/
#include <cassert>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cctype>
#include <algorithm>
using namespace std;

#include "console.h"
#include "cmdParser.h"
#include "cmdExec.h"
#include "auxil.h"
using namespace console;

extern volatile char KeepRunning;

CmdParser::CmdParser(const std::string& p) : _prompt(p), _readBufCur(0), _historyIdx(0), _dofile(nullptr){}
CmdParser::~CmdParser () {}

bool
CmdParser::regCmd(const string& cmd, const size_t nCmp, CmdExec*const e) {
    assert(e);
    // Make sure cmd hasn't been registered and won't cause ambiguity
    if (cmd.empty() || nCmp == 0 ||nCmp > cmd.size()) return false;
    for(size_t len = cmd.size(); len >= nCmp ; len --)
        if (getCmd(cmd.substr(0,len))) return false;

    string mand(nCmp,0),opt(cmd.size()-nCmp,0);
    transform(cmd.begin(),      cmd.begin()+nCmp,   mand.begin(),   ::toupper);
    transform(cmd.begin()+nCmp, cmd.end(),          opt.begin(),    ::tolower);
    e->setOptCmd(opt);
    return (_cmdMap.insert(CmdRegPair(cmd.substr(0,nCmp), e))).second;
}

//
// Parse the command from _history.back();
// Let string str = _history.back();
//
// 1. Read the command string (may contain multiple words) from the leading
//     part of str (i.e. the first word) and retrive the corresponding
//     CmdExec* from _cmdMap
//     ==> If command not found, print to cerr the following message:
//          Illegal command!! "(string cmdName)"
//     ==> return it at the end.
// 2. Call getCmd(cmd) to retrieve command from _cmdMap.
//     "cmd" is the first word of "str".
// 3. Get the command options from the trailing part of str (i.e. second
//     words and beyond) and store them in "option"
//
CmdExec*
CmdParser::parseCmd(string& option) {
    assert(_tempCmdStored.empty() && !_history.empty());
    const string& str = _history.back();

    string cmd;
    const size_t delPos = auxil::extractTok(str,cmd);
    CmdExec* cmdPtr = getCmd(cmd);
    if(!cmdPtr){
        cerr<<"Illegal command!! ("<<cmd<<")"<<endl;
        return nullptr;
    }

    try {
        option = str.substr(delPos); // might throw an out_of_range
        auxil::trimWhiteSpace(option);
    } catch (const out_of_range& e){
        option.clear();
    } return cmdPtr;
}

// cmd is a copy of the original input
///
// return the corresponding CmdExec* if "cmd" matches any command in _cmdMap
// return 0 if not found.
//
// Please note:
// ------------
// 1. The mandatory part of the command string (stored in _cmdMap) must match
// 2. The optional part can be partially omitted.
//     ==> Checked by the CmdExec::checkOptCmd(const string&) function
// 3. All string comparison are "case-insensitive".
//
CmdExec*
CmdParser::getCmd(const string& cmd) {
    for (CmdMap::const_iterator it = _cmdMap.begin();it!=_cmdMap.end();it++){
        const string mand = it->first; // mandatory part of the command
        //if (cmd.size() < mand.size()) continue;
        if (auxil::strcmpi(mand + it->second->getOptCmd(),cmd,mand.size()) ==0 )
            return it->second;
    } return nullptr;
}

// Return false on "quit" or if exception happens
CmdExecStatus
CmdParser::execOneCmd() {
    // execute the command
    if (readCmd(_dofile?*_dofile:cin)) {
        string option;
        CmdExec*const e = parseCmd(option);
        if (e) return e->exec(option);
    }
    if (KeepRunning)
        return CMD_EXEC_NOP;
    else
        return CMD_EXEC_QUIT;
}

// For each CmdExec* in _cmdMap, call its "help()" to print out the help msg.
// Print an endl at the end.
void
CmdParser::printHelps() const {
    for (CmdMap::const_iterator it = _cmdMap.begin();it!=_cmdMap.end();it++)
        it->second->help();
    cout<<endl;
}

// This function is called by pressing 'Tab'.
// It is to list the partially matched commands.
// "str" is the partial string before current cursor position. It can be
// a null string, or begin with ' '. The beginning ' ' will be ignored.
//
// Several possibilities after pressing 'Tab'
// (Let $ be the cursor position)
// 1. [Before] Null cmd
//     cmd> $
//     -----------
//     [Before] Cmd with ' ' only
//     cmd>      $
//     [After Tab]
//     ==> List all the commands, each command is printed out by:
//              cout << setw(12) << left << cmd;
//     ==> Print a new line for every 5 commands
//     ==> After printing, re-print the prompt and place the cursor back to
//          original location (including ' ')
// DONE
// 2. [Before] partially matched (multiple matches)
//     cmd> h$                         // partially matched
//     [After Tab]
//     HELp          HIStory         // List all the parially matched commands
//     cmd> h$                         // and then re-print the partial command
//     -----------
//     [Before] partially matched (multiple matches)
//     cmd> h$aaa                     // partially matched with trailing characters
//     [After Tab]
//     HELp          HIStory         // List all the parially matched commands
//     cmd> h$aaa                     // and then re-print the partial command
// DONE
// 3. [Before] partially matched (single match)
//     cmd> he$
//     [After Tab]
//     cmd> heLp $
//     -----------
//     [Before] partially matched with trailing characters (single match)
//     cmd> he$hahah
//     [After Tab]
//     cmd> heLp $hahaha
//     ==> Automatically complete on the same line
//     ==> The auto-expanded part follow the strings stored in cmd map and
//          cmd->_optCmd. Insert a space after "heLp"
//
// 4. [Before] No match
//     cmd> hek$
//     [After Tab]
//     ==> Beep and stay in the same location
//
// 5. [Before] Already matched
//     cmd> help asd$fg
//     [After] Print out the usage for the already matched command
//     Usage: HELp [(string cmd)]
//     cmd> help asd$fg DONE
//
// 6. [Before] Cursor NOT on the first word and NOT matched command
//     cmd> he haha$kk
//     [After Tab]
//     ==> Beep and stay in the same location DONE
//
void
CmdParser::listCmd(const string& str) {
    if(str.find_first_not_of(auxil::WHITESPACE)==string::npos) { // Show all commands.
        size_t i = 0;
        for(CmdMap::iterator it = _cmdMap.begin();it!=_cmdMap.end();it++,i++){
            if(i%5==0) cout<<endl;
            cout<<setw(12)<<left<<(it->first+it->second->getOptCmd());
        } return reprintCmd();
    }

    string cmd;
    if(auxil::extractTok(str,cmd)<str.size()){ // Space already inserted: search for usage.
        const CmdExec*const e = getCmd(cmd);
        if (!e) return;
        e->usage(cout << endl);
        return reprintCmd();
    }

    vector<string> match;
    for(CmdMap::iterator it = _cmdMap.begin();it!=_cmdMap.end();it++){
        const string mand = it->first + it->second->getOptCmd();
        if(auxil::strcmpi( mand, cmd, cmd.size() )==0)
            match.push_back(mand);
    }

    switch(match.size()){
        case 1:{
            for(size_t i = cmd.size(); i<match[0].size(); i++) insertChar(match[0][i]);
            insertChar(' ');
        } case 0: break;
        default:{
            for(size_t i = 0; i < match.size() ;i++){
                if(i % 5 == 0) cout<<endl;
                cout<<setw(12)<<left<<match[i];
            } return reprintCmd();
        }
    } return ;//mybeep();
}


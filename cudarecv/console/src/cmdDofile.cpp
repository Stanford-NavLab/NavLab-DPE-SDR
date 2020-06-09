
#include <cassert>
#include <fstream>
#include <stack>
#include "console.h"
#include "cmdParser.h"
//using namespace std;
using namespace console;

#define FSTACK_UB 1024 // maximum size of the file stack.

bool
CmdParser::openDofile(const std::string& dof) {
    if(_dofileStack.size()> FSTACK_UB) return false;
    _dofileStack.push(_dofile = new std::ifstream(dof.c_str()));
    if(!(*_dofile)){
        closeDofile();
        return false;
    } return _dofile;
}

// Must make sure _dofile != 0
void
CmdParser::closeDofile() {
    assert(_dofile);
    _dofile->close();
    delete _dofile;
    _dofileStack.pop();
    _dofile = _dofileStack.empty()?nullptr:_dofileStack.top();
    return;
}

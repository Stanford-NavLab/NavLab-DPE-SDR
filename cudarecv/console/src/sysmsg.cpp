
#include "console.h"
#include "cmdParser.h"

extern console::CmdParser* cmdMgr;

template <typename T> void
console::CmdParser::sysMsg (std::ostream& ostr,const T& v) {
    std::cout
        << std::string(_readBufCur + _prompt.size(),'\b')
        << std::string (_readBuf.size() + _prompt.size(),' ')
        << std::string (_readBuf.size() + _prompt.size(),'\b');
    ostr << v;

    const size_t tmp = _readBufCur;
    std::cout << _prompt << _readBuf;
    _readBufCur = _readBuf.size();
    moveBufCur(tmp);
}

template <typename T> console::os&
console::os::operator << (const T& v){
    if (!cmdMgr) _out << v;
    else cmdMgr->sysMsg(_out,v);
    return *this;
}

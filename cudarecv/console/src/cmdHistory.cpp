
#include <cassert>
//#include <string>
//#include <vector>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include "auxil.h"
#include "console.h"
#include "cmdParser.h"
//using namespace std;
using namespace console;


// This function adds the string in _readBuf to the _history.
// The size of _history may or may not change. Depending on whether
// there is a temp history string.
//
// 1. Remove ' ' at the beginning and end of _readBuf
// 2. If not a null string, add string to _history.
//     Be sure you are adding to the right entry of _history.
// 3. If it is a null string, don't add anything to _history.
// 4. Make sure to clean up "temp recorded string" (added earlier by up/pgUp,
//     and reset _tempCmdStored to false
// 5. Reset _historyIdx to _history.size() // for future insertion
// 6. Reset _readBufPtr and _readBufEnd to _readBuf
// 7. Make sure *_readBufEnd = 0 ==> _readBuf becomes null string
//
bool
CmdParser::addHistory() {
    _tempCmdStored.clear();

    if(!auxil::trimWhiteSpace(_readBuf)) return false;
    _history.push_back(_readBuf);
    _historyIdx = _history.size();
    return true;
}


// 1. Replace current line with _history[_historyIdx] on the screen
// 2. Set _readBufPtr and _readBufEnd to end of line
//
// [Note] Do not change _history.size().
//
void
CmdParser::retrieveHistory() {
    deleteLine();
    _readBuf = _historyIdx < _history.size()? _history[_historyIdx] :_tempCmdStored;
    std::cout << _readBuf;
    _readBufCur = _readBuf.size();
}

void
CmdParser::moveToHistory(ptrdiff_t index) {
    if (_history.empty()) return;

    // Move back to current.
    if (index >= ptrdiff_t(_history.size()))
        _historyIdx = _history.size();
    else if (_historyIdx == _history.size()){// Leaving current
        assert(_tempCmdStored.empty()); // TODO: change to throw
        _tempCmdStored = _readBuf;
    } _historyIdx = (index < 0 ? 0 : index) ;
    retrieveHistory();
}

void
CmdParser::printHistory(const size_t nPrint) const {
    assert(_tempCmdStored.empty()); // TODO: change to throw
    if (_history.empty()) {
        std::cout << "Empty command history!!" << std::endl;
        return;
    }

    const size_t sz = _history.size();
    for (size_t i = nPrint?(sz - std::min(nPrint,sz)):0; i < sz; i++)
        std::cout << std::setw(5) << std::right<< i << ": " << std::left << _history[i] << std::endl;
}


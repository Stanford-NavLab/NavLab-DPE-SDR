/****************************************************************************
  FileName     [ cmdParser.cpp ]
  PackageName  [ cmd ]
  Synopsis     [ Define command line reader member functions ]
  Author       [ Chung-Yang (Ric) Huang ]
  Copyright    [ Copyleft(c) 2007-2012 LaDs(III), GIEE, NTU, Taiwan ]
****************************************************************************/
#include <cassert>
#include <string>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include "console.h"
#include "cmdParser.h"
using namespace std;
using namespace console;

extern volatile char KeepRunning;

bool
CmdParser::moveBufCur(const size_t pos){
    if (pos > _readBuf.size()) throw out_of_range("pos");
    for(;pos<_readBufCur;_readBufCur--)
        cout<<'\b';
    for(;pos>_readBufCur;_readBufCur++)
        cout<<_readBuf[_readBufCur];
    return true;
}

bool
CmdParser::deleteChar() {
    if(_readBufCur==_readBuf.size()) return false;
    const size_t prev = _readBufCur;
    cout << _readBuf.erase(_readBufCur,1).substr(_readBufCur) << ' ';
    _readBufCur = _readBuf.size() + 1;
    moveBufCur(prev);
    return true;
}

void
CmdParser::insertChar(const char ch, const size_t rep) {
    const size_t prev = _readBufCur;
    cout << _readBuf.insert(_readBufCur,rep,ch).substr(_readBufCur);
    _readBufCur = _readBuf.size();
    moveBufCur(prev + rep);
    return;
}

void
CmdParser::deleteLine() {
    moveBufCur(0);
    cout << string(_readBuf.size(),' ');
    _readBufCur = _readBuf.size();
    _readBuf.clear();
    moveBufCur(0);
    return;
}

void
CmdParser::reprintCmd() {
     const size_t tmp = _readBufCur;
     cout << endl << _prompt << _readBuf;
     _readBufCur = _readBuf.size();
     moveBufCur(tmp);
}

CmdParser::ParseChar
CmdParser::getChar(istream& istr) const {
    char ch = getKey(istr);

    if (istr.eof()) return ParseChar(INPUT_END_KEY);

    switch (ch) {
        // Simple keys: one code for one key press
        // -- The following should be platform-independent
        case LINE_BEGIN_KEY:  // Ctrl-a
        case LINE_END_KEY:    // Ctrl-e
        case INPUT_END_KEY:   // Ctrl-d
        case TAB_KEY:         // tab('\t') or Ctrl-i
        case NEWLINE_KEY:     // enter('\n') or ctrl-m
        case BACK_SPACE_KEY:
        case BACK_SPACE_CHAR:
            return ParseChar(ch);

        // keyboard mapping for special keys.
        // -- The following simple/combo keys are platform-dependent
        //     You should test to check the returned codes of these key presses
        // -- You should either modify the "enum ParseChar" definitions in
        //     "charDef.h", or revise the control flow of the "case ESC" below
        //
        // -- You need to handle:
        //     { BACK_SPACE_KEY, ARROW_UP/DOWN/RIGHT/LEFT,
        //        HOME/END/PG_UP/PG_DOWN/INSERT/DELETE }
        //
        // Combo keys: multiple codes for one key press
        // -- Usually starts with ESC key, so we check the "case ESC"
        // case ESC_KEY:

        // For the remaining printable and undefined keys

        case ESC_KEY: {
            char combo = getKey(istr);
            // Note: ARROW_KEY_INT == MOD_KEY_INT, so we only check MOD_KEY_INT
            if (combo == char(MOD_KEY_INT) || combo == char (HE_KEY_INT)) {
                char key = getKey(istr);
                if ((key >= char(MOD_KEY_BEGIN)) && (key <= char(MOD_KEY_END)) && (key != char(FALSY_END))) {
                    if (getKey(istr) == MOD_KEY_DUMMY)
                        return ParseChar(int(key) + MOD_KEY_FLAG);
                    else return ParseChar(UNDEFINED_KEY);

                } else if ((key >= char(ARROW_KEY_BEGIN)) && (key <= char(ARROW_KEY_END)))
                    return ParseChar(int(key) + ARROW_KEY_FLAG);

                else if(key == char(HOME_KEY)||key == char(END_KEY))
                    return ParseChar (int(key) + MOD_KEY_FLAG);

                else return ParseChar(UNDEFINED_KEY);
            } /*else if (combo == char(HE_KEY_INT)){
                char key = getKey(istr);
                if(key == char(HOME_KEY)||key == char(END_KEY))
                    return ParseChar (int(key) + MOD_KEY_FLAG);
                else return ParseChar(UNDEFINED_KEY);
            } */else { /*mybeep();*/ return getChar(istr); }
        }
        default:
            if (isprint(ch)) return ParseChar(ch);
            else return ParseChar(UNDEFINED_KEY);
    } return ParseChar(UNDEFINED_KEY);
}

bool
CmdParser::readCmd(istream& istr) {
//    resetBuf();
//    printPrompt();
    _readBuf.clear();
    _readBufCur = 0;
    cout << _prompt;

    while (KeepRunning) {
        ParseChar pch = getChar(istr);
        if (pch == INPUT_END_KEY && _dofile) {
           closeDofile();
           return false;
        }
        switch(pch) {
            case LINE_BEGIN_KEY :
            case HOME_KEY       : moveBufCur(0); break;
            case LINE_END_KEY   :
            case END_KEY        : moveBufCur(_readBuf.size()); break;
            case BACK_SPACE_KEY :
            case BACK_SPACE_CHAR:{
                if (_readBufCur) moveBufCur(_readBufCur - 1);
                else break;
            }
            case DELETE_KEY     : deleteChar(); break;
            case NEWLINE_KEY    : {
                cout << char(NEWLINE_KEY);
                return addHistory();
            }
            case ARROW_UP_KEY   : moveToHistory(_historyIdx - 1); break;
            case ARROW_DOWN_KEY : moveToHistory(_historyIdx + 1); break;
            case ARROW_RIGHT_KEY: moveBufCur(_readBufCur <= _readBuf.size()?_readBufCur+1:_readBuf.size()); break;
            case ARROW_LEFT_KEY : moveBufCur(_readBufCur > 0 ? _readBufCur - 1:0); break;
            case PG_UP_KEY      : moveToHistory(_historyIdx - PG_OFFSET); break;
            case PG_DOWN_KEY    : moveToHistory(_historyIdx + PG_OFFSET); break;
            case TAB_KEY        : listCmd(_readBuf); break;//insertChar(' ',TAB_POSITION-(_readBufCur%TAB_POSITION)); break;
            case INSERT_KEY     : // not yet supported; fall through to UNDEFINE
            case UNDEFINED_KEY  : //mybeep();
            break;
            default: insertChar(char(pch)); break;
        }
    } return false;/// newline;
}


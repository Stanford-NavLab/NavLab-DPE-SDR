
#ifndef _CMD_READER_H
#define _CMD_READER_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stack>
#include <cstddef>
#include <map>
#include <climits>
#include "console.h"
#include "cmdExec.h"

#define TAB_POSITION 8
#define PG_OFFSET    10


//----------------------------------------------------------------------
//     Base class : Console
//----------------------------------------------------------------------

class console::CmdParser {

public:
    enum ParseChar:int;

    /* cmdParser.cpp */
    CmdParser(const std::string& p);
    virtual ~CmdParser ();

    bool          regCmd         (const std::string&,const size_t nCmp, CmdExec*const e);
    CmdExec*      getCmd         (const std::string&);
    void          printHelps     () const;
    CmdExecStatus execOneCmd();

    /* cmdDofile.cpp */
    bool    openDofile     (const std::string& dof);
    void    closeDofile    ();

    void    printHistory   (const size_t nPrint = 0) const;
    template <typename T> void
            sysMsg         (std::ostream&,const T&);

private:
    typedef std::map  <const std::string, console::CmdExec*>  CmdMap;
    typedef std::pair <const std::string, console::CmdExec*>  CmdRegPair;

    /* cmdParser.cpp */
    CmdExec* parseCmd   (std::string& option);
    void     listCmd    (const std::string& cmd);

    // Private member functions
    //void resetBuf       ();
    //void printPrompt    () const;

    // Helper functions
    /* cmdReader.cpp */
    bool moveBufCur     (const size_t);
    bool deleteChar     ();
    void insertChar     (const char ch, const size_t rep = 1);
    void deleteLine     ();
    void reprintCmd     ();
    bool      readCmd   (std::istream&);
    ParseChar getChar   (std::istream&) const;

    /* cmdHistory.cpp */
    void moveToHistory  (ptrdiff_t index);
    bool addHistory     ();
    void retrieveHistory();


    // Data members
    const std::string           _prompt;
    size_t                      _readBufCur,_historyIdx;
    std::string                 _readBuf,_tempCmdStored;
    std::ifstream*              _dofile;
    std::vector<std::string>    _history;
    std::stack<std::ifstream*>  _dofileStack;
    CmdMap                      _cmdMap;                // std::map from std::string to command

public:
    enum ParseChar :int{
        LINE_BEGIN_KEY   = 1,          // ctrl-a
        LINE_END_KEY     = 5,          // ctrl-e
        INPUT_END_KEY    = 4,          // ctrl-d
        TAB_KEY          = int('\t'),  // tab('\t') or Ctrl-i
        NEWLINE_KEY      = int('\n'),  // enter('\n') or ctrl-m
        ESC_KEY          = 27,         // Not printable; used for combo keys

        // -- The following simple/combo keys are platform-dependent
        //     You should test to check the returned codes of these key presses
        // -- Use "testAsc.cpp" to test
        //
        // [FLAG bit for combo keys]
        // -- Added to the returned ParseChar of combo keys
        // -- i.e. The returned ParseChar will be "ComboKeyEnum + FLAG bit"
        // -- This is to avoid the collision with the ASCII codes of regular keys
        // -- Feel free to add/remove/modify on your own
        //
        // [Intermediate keys for combo keys]
        // -- Intermediate keys are the common parts of combo keys
        //
        BACK_SPACE_KEY      = 127 ,

        //
        // -- Arrow keys: 27 -> 91 -> {UP=65, DOWN=66, RIGHT=67, LEFT=68}
        ARROW_KEY_FLAG   = 1 << 8,
        ARROW_KEY_INT    = 91,
        ARROW_UP_KEY     = 65 + ARROW_KEY_FLAG,
        ARROW_DOWN_KEY   = 66 + ARROW_KEY_FLAG,
        ARROW_RIGHT_KEY  = 67 + ARROW_KEY_FLAG,
        ARROW_LEFT_KEY   = 68 + ARROW_KEY_FLAG,
        ARROW_KEY_BEGIN  = ARROW_UP_KEY,
        ARROW_KEY_END    = ARROW_LEFT_KEY,

        //
        // -- MOD keys: 27 -> 91 -> {49-54} -> 126
        //     MOD_KEY = { INSERT, DELETE, HOME, END, PgUp, PgDown }
        //

        MOD_KEY_FLAG      = 1 << 9,
        MOD_KEY_INT       = 91,
        INSERT_KEY        = 50 + MOD_KEY_FLAG,
        DELETE_KEY        = 51 + MOD_KEY_FLAG,
        PG_UP_KEY         = 53 + MOD_KEY_FLAG,
        PG_DOWN_KEY       = 54 + MOD_KEY_FLAG,
        MOD_KEY_BEGIN     = INSERT_KEY,
        MOD_KEY_END       = PG_DOWN_KEY,
        MOD_KEY_DUMMY     = 126,

        HE_KEY_INT        = 79,
        HOME_KEY          = 72 + MOD_KEY_FLAG,
        END_KEY           = 70 + MOD_KEY_FLAG,
        FALSY_END         = 52 + MOD_KEY_FLAG,
        //
        // [For undefined keys]
        UNDEFINED_KEY     = INT_MAX,

        // For output only, you don't need to modify this part
        BEEP_CHAR         = 7,
        BACK_SPACE_CHAR   = 8,

        // dummy end
        PARSE_CHAR_END
    };

private:

};

#endif // CMD_PARSER_H

/****************************************************************************
  FileName     [ cmdParser.h ]
  PackageName  [ cmd ]
  Synopsis     [ Define class Console ]
  Author       [ Chung-Yang (Ric) Huang ]
  Copyright    [ Copyleft(c) 2007-2013 LaDs(III), GIEE, NTU, Taiwan ]
****************************************************************************/

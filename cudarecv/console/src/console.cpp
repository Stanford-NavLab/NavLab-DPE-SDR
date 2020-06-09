
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cctype>
#include <termios.h>
using namespace std;

#include "console.h"
#include "cmdParser.h"
//#include "cmdCommon.h"
using namespace console;

static struct termios stored_settings;

static void set_keypress(void) {
    struct termios new_settings;
    tcgetattr(0,&stored_settings);
    new_settings = stored_settings;
    new_settings.c_lflag &= (~ICANON);
    new_settings.c_lflag &= (~ECHO);
    new_settings.c_cc[VTIME] = 0;
    tcgetattr(0,&stored_settings);
    new_settings.c_cc[VMIN] = 1;
    tcsetattr(0,TCSANOW,&new_settings);
}

char
console::getKey(istream& istr) {
    char ch;
    set_keypress();
    istr.unsetf(ios_base::skipws);
    istr >> ch;
    istr.setf(ios_base::skipws);
    tcsetattr(0,TCSANOW,&stored_settings);
    //reset_keypress();
    return ch;
}


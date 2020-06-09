
#include <iostream>
#include <getopt.h>
#include <cstring>
#include <stdint.h>
#include <cstdlib>
#include "startoptions.h"
#include "auxil.h"

using namespace std;

#define OPT_STR    "f:"

void StartOptions::Usage(char* program) {
    cout << "Usage: " << program << "[-i] [-f filename]" << endl;
    cout << "       " << "--no-console" << endl;
    cout << "       " << "--filename filename" << endl;
    cout << "       " << "--samplerate 2000000" << endl;
    cout << "       " << "--carrierfrequency 0" << endl;
    cout << "       " << "--skip 0" << endl;
    cout << "       " << "--fileformat 0" << endl;
}

void StartOptions::ParseOptions(int argc, char* argv[]) {
    int c;

    while (1) {
        static struct option longOptions[] = {
            {"no-console", no_argument, 0, '5'},
            {"filename", required_argument, 0, 'f'},
            {"samplerate", required_argument, 0, '1'},
            {"carrierfrequency", required_argument, 0, '2'},
            {"skip", required_argument, 0, '3'},
            {"fileformat", required_argument, 0, '4'},
            {0, 0, 0, 0}
        };

        int optionIndex = 0;

        c = getopt_long(argc, argv, OPT_STR, longOptions, &optionIndex);

        if (c == -1) break;    // end of options

        switch (c) {
        case '5':
            console = false;
            break;
        case 'f':
            fromFile = true;
            filename = optarg;
            //cout << "from file: " << filename << endl;
            break;
        case '1':
            //parse sample rate
            sampleRate = ParseFreq(optarg);
            //cout << "sample rate: " << sampleRate << endl;
            break;
        case '2':
            carrierFreq = ParseFreq(optarg);
            //cout << "carrier frequency: " << carrierFreq << endl;
            //parse carrier frequency
            break;
        case '3':
            //cout << "skip: " << optarg << endl;
            //parse skip
            break;
        case '4':
            //cout << "file format: " << optarg << endl;
            //parse fileformat
            break;
        case '?':
            /* getopt_long already printed an error message. */
            error = true;
            Usage(argv[0]);
            return;
            break;
        default:
            cerr << "Error: option " << c
                 << " not recognized. Argument: " << optarg << endl;
            error = true;
            Usage(argv[0]);
            return;
        }
    }

    //cout << "console: " << console << endl;

}

StartOptions::StartOptions(int argc, char* argv[]) {
    InitDefaults();
    ParseOptions(argc, argv);
}

void StartOptions::InitDefaults() {
    error = false;
    console = true;
    fromFile = false;
    sampleRate = 2000000;
    carrierFreq = 0;
    skip = 0;
    fileFormat = 0;
    filename.clear();
}

int StartOptions::ParseFreq(char* arg) {
    double val = 0.0;
    string units;
    SplitArgument(arg, val, units);
    if (units.length() > 0){
        if (((auxil::strcmpi(units, string("MHz")) == 0) && (units.at(0) == 'M')) ||
                (units.compare("M") == 0))
            val *= 1000000.0;
        else if ((auxil::strcmpi(units, string("GHz")) == 0) ||
                (auxil::strcmpi(units, string("G")) == 0))
            val *= 1000000000.0;
        else if ((auxil::strcmpi(units, string("kHz")) == 0) ||
                (auxil::strcmpi(units, string("k")) == 0))
            val *= 1000.0;
        else
            return -1;
    }
    val += 0.5;
    return (int)val;
}

int StartOptions::SplitArgument(char* arg, double &val, string &units) {
    string str = arg;
    auxil::trimWhiteSpace(str);
    size_t split = str.find_first_not_of("1234567890.");
    string num;
    if (split >= str.length()){
        units.clear();
        num = str;
    } else {
        num = str.substr(0,split);
        str = str.substr(split, string::npos);
        auxil::trimWhiteSpace(str);
        units = str;
    }
    val = atof(num.c_str());
    return 0;
}


#ifndef INC__STARTOPTIONS_H_
#define INC__STARTOPTIONS_H_

#include <string>

//using namespace std;

class StartOptions {

    public:

        StartOptions(int argc, char* argv[]);

        bool console;                   /**< console mode */
        bool fromFile;                  /**< samples from file */
        unsigned int sampleRate;        /**< in Hz */
        unsigned int carrierFreq;       /**< in Hz */
        unsigned long skip;             /**< in number of samples */
        unsigned char fileFormat;
        std::string filename;

        bool error;

    private:

        /** \brief Initializes all options to defaults. */
        void InitDefaults(void);

        /** \brief Parses all of argv[] into options. */
        void ParseOptions(int argc, char* argv[]);

        /** \brief Prints usage instructions. */
        void Usage(char* program);

        /** \brief Parses frequency string into integer.
         *
         *  Parses possible strings such as 12.3GHz into integer of value Hz.
         *  Possible units are GHz, MHz, kHz.
         *  return frequency in Hz. -1 on error.
         */
        int ParseFreq(char* arg);

        /** \brief Splits an input string into its value and units. */
        int SplitArgument(char* arg, double &val, std::string &units);

};

#endif

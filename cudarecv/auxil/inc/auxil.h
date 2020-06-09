
#ifndef INC__AUXIL_H_
#define INC__AUXIL_H_

#include <iostream>
#include <string>
#include <pthread.h>
#include <map>

/** \brief  Return the positive modulo (wrap around) of a % b (like Python) */
#define POSMOD(a,b) ((((a)%(b))+(b))%(b))

/**
 * I wanted to have a collection of helper functions but dont want to saturate
 * the global namespace.
 */
namespace auxil {
    const char* const WHITESPACE = " \t\n\v\f\r";
    /** \brief Case-insensitive version of strcmp. */
    int strcmpi(std::string templt, std::string target,const size_t reqN = 0) ;

    /** \brief Function name says it all. */
    size_t trimWhiteSpace(std::string &str);
    size_t extractTok(const std::string& str, std::string& tok, const size_t pos = 0, const char* del = WHITESPACE) ;

    uint64_t ClockGetTime();
    unsigned int roundUpToNextPowerOfTwo(unsigned int x);
    void   create_pthread_attr(pthread_attr_t&,const int);
	double str2num(const char *s, int i, int n);


    // define a data structure to store time and memory usage
    struct TmStat {
	    long vmSize; // in kilobytes
	    long vmPeak; // in kilobytes
	    long vmDiff; // in kilobytes
	    long rTime;  // real time in micro seconds
	    long uTime;  // user time in micro seconds
	    long sTime;  // system time in micro seconds
    };

    class TmUsage {
	    public:
		    TmUsage();
		    ~TmUsage();

	    bool totalStart();     // start the total timer.  this is called only once.
	    bool periodStart();    // start the period timer. this can be called many times
	    bool getTotalUsage(TmStat &st) const;   // get the total tm usage
	    bool getPeriodUsage(TmStat &st) const;  // get the priod tm usage
    void printUsage(std::ostream& os,const TmStat&st) const;

	    private:
	    TmStat tStart_;
	    TmStat pStart_;
	    bool checkUsage(TmStat &st) const;      //
    };
    
    
    // From: https://stackoverflow.com/questions/5009071/c-templated-function-error-while-trying-to-throw-exception
    class Invalid_State_Exception : public std::runtime_error
    {
    public:
        Invalid_State_Exception(const std::string& msg) :
            std::runtime_error(msg) { }
    };
    
    // From: https://stackoverflow.com/questions/726664/string-to-enum-in-c
    template <typename T>
    class EnumParser
    {
	    std::map<std::string, T> enumMap;
	
	public:
	    EnumParser(){};

	    T ParseSomeEnum(const std::string &value);
    };
	
    #include "auxil.tpp"

}

#endif

// **************************************************************************
// File       [ tm_usage.h ]
// Author     [ littleshamoo ]
// Synopsis   [ Get CPU time and Memory usage ]
// Date       [ Ver 2.0 started 2010/03/23 ]
// History    [ created TmStat structure to store usage ]
// **************************************************************************



#include <algorithm>
#include <string.h>
#include "auxil.h"
#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif
using namespace std;

//#define

int auxil::strcmpi(string templt, string target,const size_t reqN) {
    std::transform(templt.begin(), templt.end(), templt.begin(), ::tolower);
    std::transform(target.begin(), target.end(), target.begin(), ::tolower);
    if (reqN == 0 || templt.size() < reqN || target.size() < reqN || target.size() > templt.size()) return templt.compare(target);
    return templt.substr(0,target.size()).compare(target);
}

size_t auxil::trimWhiteSpace(string &str) {
//    if (str.length() == 0) return;
 //   size_t first = str.find_first_not_of(' ');
//    size_t last = str.find_last_not_of(' ');
    const size_t pos = str.find_first_not_of(WHITESPACE);
    const size_t len = str.find_last_not_of (WHITESPACE) - pos + 1;
    if (pos == string::npos){
        str.clear();
        return 0;
    } return (str = str.substr(pos,len)).size();
}

size_t auxil::extractTok(const string& str, string& tok, const size_t pos, const char* del){
    size_t begin = str.find_first_not_of(del, pos);
    if (begin == string::npos) { tok.clear(); return begin; }
    size_t end = str.find_first_of(del, begin);
    tok = str.substr(begin, end - begin);
    return end;
}

void auxil::create_pthread_attr (pthread_attr_t& attr,const int prio){
    pthread_attr_init (&attr);
    pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_JOINABLE);

    struct sched_param param;
    int minPrio = sched_get_priority_min(SCHED_RR);
    int maxPrio = sched_get_priority_max(SCHED_RR);

    param.sched_priority = (maxPrio - minPrio) *
                           prio / 32 + minPrio;

    pthread_attr_setschedpolicy(&attr, SCHED_RR);
    pthread_attr_setschedparam (&attr, &param);
}


/* string to number ------------------------------------------------------------
* convert substring in string to number
* args   : char   *s        I   string ("... nnn.nnn ...")
*          int    i,n       I   substring position and width
* return : converted number (0.0:error)
*
* Taken from rtkcmn.c in rtklib
*-----------------------------------------------------------------------------*/
double auxil::str2num(const char *s, int i, int n)
{
    double value;
    char str[256],*p=str;

    if (i<0||(int)strlen(s)<i||(int)sizeof(str)-1<n) return 0.0;
    for (s+=i;*s&&--n>=0;s++) *p++=*s=='d'||*s=='D'?'E':*s; *p='\0';
    return sscanf(str,"%lf",&value)==1?value:0.0;
}



#include <time.h>
#include <sys/timeb.h>
// needs -lrt (real-time lib)
// 1970-01-01 epoch UTC time, 1 mcs resolution (divide by 1M to get time_t)
uint64_t auxil::ClockGetTime() {
    timespec ts;
#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    ts.tv_sec = mts.tv_sec;
    ts.tv_nsec = mts.tv_nsec;
#else
    clock_gettime(CLOCK_REALTIME, &ts);
#endif
    return (uint64_t)ts.tv_sec * 1000000LL + (uint64_t)ts.tv_nsec / 1000LL;
} // http://stackoverflow.com/questions/275004



unsigned int auxil::roundUpToNextPowerOfTwo(unsigned int x)
{
    x--;
    x |= x >> 1;  // handle  2 bit numbers
    x |= x >> 2;  // handle  4 bit numbers
    x |= x >> 4;  // handle  8 bit numbers
    x |= x >> 8;  // handle 16 bit numbers
    x |= x >> 16; // handle 32 bit numbers
    x++;
    
    return x;
} // https://bits.stephan-brumme.com/roundUpToNextPowerOfTwo.html





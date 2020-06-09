
#include "converters.h"


/* add time --------------------------------------------------------------------
* add time to gtime_t struct
* args   : gtime_t t        I   gtime_t struct
*          double sec       I   time to add (s)
* return : gtime_t struct (t+sec)
*
* Taken from rtkcmn.c
*-----------------------------------------------------------------------------*/
dsp::utils::gtime_t dsp::utils::timeadd(dsp::utils::gtime_t t, double sec)
{
    double tt;

    t.sec+=sec; tt=floor(t.sec); t.time+=(int)tt; t.sec-=tt;
    return t;
}

/* time difference -------------------------------------------------------------
* difference between gtime_t structs
* args   : gtime_t t1,t2    I   gtime_t structs
* return : time difference (t1-t2) (s)
*
* Taken from rtkcmn.c
*-----------------------------------------------------------------------------*/
double dsp::utils::timediff(dsp::utils::gtime_t t1, dsp::utils::gtime_t t2)
{
    return difftime(t1.time,t2.time)+t1.sec-t2.sec;
}

/* string to time --------------------------------------------------------------
* convert substring in string to gtime_t struct
* args   : char   *s        I   string ("... yyyy mm dd hh mm ss ...")
*          int    i,n       I   substring position and width
*          gtime_t *t       O   gtime_t struct
* return : status (0:ok,0>:error)
*
* Taken from rtkcmn.c in rtklib
*-----------------------------------------------------------------------------*/
int dsp::utils::str2time(const char *s, int i, int n, dsp::utils::gtime_t *t)
{
    double ep[6];
    char str[256],*p=str;

    if (i<0||(int)strlen(s)<i||(int)sizeof(str)-1<i) return -1;
    for (s+=i;*s&&--n>=0;) *p++=*s++; *p='\0';
    if (sscanf(str,"%lf %lf %lf %lf %lf %lf",ep,ep+1,ep+2,ep+3,ep+4,ep+5)<6)
        return -1;
    if (ep[0]<100.0) ep[0]+=ep[0]<80.0?2000.0:1900.0;
    *t=dsp::utils::epoch2time(ep);
    return 0;
}


/* adjust time considering week handover
 *
 * Taken from rinex.c
 */
dsp::utils::gtime_t dsp::utils::adjweek(dsp::utils::gtime_t t, dsp::utils::gtime_t t0)
{
    double tt=timediff(t,t0);
    if (tt<-302400.0) return dsp::utils::timeadd(t, 604800.0);
    if (tt> 302400.0) return dsp::utils::timeadd(t,-604800.0);
    return t;
}

/* convert calendar day/time to time -------------------------------------------
* convert calendar day/time to gtime_t struct
* args   : double *ep       I   day/time {year,month,day,hour,min,sec}
* return : gtime_t struct
* notes  : proper in 1970-2037 or 1970-2099 (64bit time_t)
*
* Taken from rtkcmn.c
*-----------------------------------------------------------------------------*/
dsp::utils::gtime_t dsp::utils::epoch2time(const double *ep)
{
    const int doy[]={1,32,60,91,121,152,182,213,244,274,305,335};
    dsp::utils::gtime_t time={0};
    int days,sec,year=(int)ep[0],mon=(int)ep[1],day=(int)ep[2];

    if (year<1970||2099<year||mon<1||12<mon) return time;

    /* leap year if year%4==0 in 1901-2099 */
    days=(year-1970)*365+(year-1969)/4+doy[mon-1]+day-2+(year%4==0&&mon>=3?1:0);
    sec=(int)floor(ep[5]);
    time.time=(time_t)days*86400+(int)ep[3]*3600+(int)ep[4]*60+sec;
    time.sec=ep[5]-sec;
    return time;
}



/* gps time to time ------------------------------------------------------------
* convert week and tow in gps time to gtime_t struct
* args   : int    week      I   week number in gps time
*          double sec       I   time of week in gps time (s)
* return : gtime_t struct
*
* Taken from rtkcmn.c
*-----------------------------------------------------------------------------*/
dsp::utils::gtime_t dsp::utils::gpst2time(int week, double sec)
{
    dsp::utils::gtime_t t=dsp::utils::epoch2time(dsp::utils::GPST0);

    if (sec<-1E9||1E9<sec) sec=0.0;
    t.time+=86400*7*week+(int)sec;
    t.sec=sec-(int)sec;
    return t;
}

/* time to gps time ------------------------------------------------------------
* convert gtime_t struct to week and tow in gps time
* args   : gtime_t t        I   gtime_t struct
*          int    *week     IO  week number in gps time (NULL: no output)
* return : time of week in gps time (s)
*
* Taken from rtkcmn.c
*-----------------------------------------------------------------------------*/
double dsp::utils::time2gpst(dsp::utils::gtime_t t, int *week)
{
    dsp::utils::gtime_t t0=dsp::utils::epoch2time(dsp::utils::GPST0);
    time_t sec=t.time-t0.time;
    int w=(int)(sec/(86400*7));

    if (week) *week=w;
    return (double)(sec-w*86400*7)+t.sec;
}

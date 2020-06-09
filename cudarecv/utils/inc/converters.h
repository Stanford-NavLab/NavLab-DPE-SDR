
#ifndef INC__CONVERTERS_H_
#define INC__CONVERTERS_H_

#include <iostream>
#include <string>
#include <pthread.h>
#include <map>
#include <ephhelper.h>
#include <math.h>
#include <string.h>


/**
 * I wanted to have a collection of helper functions but dont want to saturate
 * the global namespace.
 */
namespace dsp { namespace utils {


	gtime_t timeadd(gtime_t t, double sec);
	double timediff(gtime_t t1, gtime_t t2);
	int str2time(const char *s, int i, int n, gtime_t *t);
	gtime_t adjweek(gtime_t t, gtime_t t0);
	gtime_t epoch2time(const double *ep);
	gtime_t gpst2time(int week, double sec);
	double time2gpst(gtime_t t, int *week);


}
}


#endif

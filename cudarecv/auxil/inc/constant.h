
#ifndef _CONSTANT_H
#define _CONSTANT_H

#include <cmath>

#define WGS84_c             (299792458.0)
#define WGS84_MU            (3.986005e14)
#define WGS84_OMEGA_DOT_E   (7.2921151467e-5)
#define WGS84_a             (6378137.0)
#define WGS84_f             (1.0/298.257223563)// TODO: verify the number
#define REL_F               (-4.442807633E-10)

#define WGS84_b             (WGS84_a * (1.0-WGS84_f))
#define WGS84_e             sqrt(1-WGS84_b*WGS84_b/WGS84_a/WGS84_a)
#define WGS84_ep            sqrt(WGS84_a*WGS84_a/WGS84_b/WGS84_b-1)

#define T_CA                (1e-3)
#define F_CA                (1.023e6)
#define L_CA				(1023)
#define F_L1                (1575.42e6)
#define WGS84_PI            (3.1415926535898)



#endif

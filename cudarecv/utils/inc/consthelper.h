#ifndef __INC_UTILS_CONSTHELPER_H_
#define __INC_UTILS_CONSTHELPER_H_


#define CONST_C 			(299792458.0)
#define CONST_PI 			(3.1415926535898)
#define CONST_2PI			(6.2831853071796)		// Save some multiplications on the GPU
#define CONST_F_L1			(1.57542e9)
#define CONST_F_L2			(1.22760e9)

#define CONST_F_CA			(1.023e6)
#define CONST_L_CA			(1023.0)
#define CONST_T_CA			(0.001)
#define CONST_PRN_MAX		(37)
#define CONST_GPS_WEEK_SIZE	(604800)

/** \brief Constants useful for GPS processing */
#define CONST_MU 			(3.986005e14) 			/**< earth's universal gravitation parameter */
#define CONST_F 			(-4.442807633e-10)      /**< earth's relativistic correction term */
#define CONST_OEDot 		(7.2921151467e-5)		/**< earth's sidereal rotation rate */

/** \brief Ellipsoid conversion constants */
#define CONST_WGS84_A 		(6378137.0)
#define CONST_WGS84_B		(6356752.314245)
#define CONST_WGS84_E		(0.08181919084262149)
#define CONST_WGS84_EP		(0.08209443794969568)
#define CONST_WGS84_INVF 	(298.257223563)

// Formulae for CONSTs -- better to use pre-calc'd values to ensure they're not being computed at runtime
//#define CONST_WGS84_B		(CONST_WGS84_A * (1.0 - (1.0/CONST_WGS84_INVF)))
//#define CONST_WGS84_E		(sqrt(CONST_WGS84_A*CONST_WGS84_A - CONST_WGS84_B*CONST_WGS84_B)/(CONST_WGS84_A*CONST_WGS84_A))
//#define CONST_WGS84_EP		(sqrt(CONST_WGS84_A*CONST_WGS84_A - CONST_WGS84_B*CONST_WGS84_B)/(CONST_WGS84_B*CONST_WGS84_B))

/*
namespace dsp { namespace utils {


	static const double C = 299792458.0;            // < speed of light in meters per second
	static const double PI = 3.1415926535898;       // < GPS pi
	static const double F_L1 = 1.57542e9;           // < GPS L1 frequency (Hz)
	static const double F_L2 = 1.22760e9;           // < GPS L2 frequency (Hz)

	static const double F_CA = 1.023e6; 			// < CA code frequency (chips/s)
	static const double L_CA = 1023.0;  			// < CA code chips in a period (chips)
	static const double T_CA = 0.001;   			// < CA code period (T_CA = L_CA/f_CA)
	static const double PRN_MAX = 37;				// < Range of PRNs (1-37) that we have codes stored for

}
}
*/

#endif

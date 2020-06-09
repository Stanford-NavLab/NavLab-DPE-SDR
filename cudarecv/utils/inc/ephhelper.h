#ifndef __INC_UTILS_EPHHELPER_H_
#define __INC_UTILS_EPHHELPER_H_

#include "consthelper.h"
#include <vector>

//#define PI              3.141592654
#define MAX_CHAN        1

namespace dsp { namespace utils {


	const static double GPST0[]={1980,1, 6,0,0,0}; 	/* gps time reference */

	static const double URA_EPH[]={         		/* ura values (ref [3] 20.3.3.3.1.1) */
	    2.4,3.4,4.85,6.85,9.65,13.65,24.0,48.0,96.0,192.0,384.0,768.0,1536.0,
	    3072.0,6144.0,0.0
	};
	// Ref [3]: https://www.gps.gov/technical/icwg/IS-GPS-800A.pdf -> See page 43 (per bottom of page)

	/** \brief RTKLIB RINEX ephemerides types */
	static const int TSYS_GPS = 0;
	static const int TSYS_UTC = 1;
	static const int TSYS_GLO = 2;
	static const int TSYS_GAL = 3;
	static const int TSYS_QZS = 4;
	static const int TSYS_CMP = 5;
	static const int TSYS_IRN = 6;

	/** \brief RTKLIB RINEX ephemeris mask types */
	static const unsigned int SYS_NONE	= 0x00;
	static const unsigned int SYS_GPS	= 0x01;
	static const unsigned int SYS_SBS	= 0x02;
	static const unsigned int SYS_GLO	= 0x04;
	static const unsigned int SYS_GAL	= 0x08;
	static const unsigned int SYS_QZS	= 0x10;
	static const unsigned int SYS_CMP	= 0x20;
	static const unsigned int SYS_IRN	= 0x40;
	static const unsigned int SYS_LEO	= 0x80;
	static const unsigned int SYS_ALL	= 0xFF;

	/** \brief RTKLIB RINEX satellite constellation constants */
	static const unsigned int MINPRNGPS = 1;
	static const unsigned int MAXPRNGPS = 32;
	static const unsigned int NSATGPS = (MAXPRNGPS-MINPRNGPS+1);
	static const unsigned int NSYSGPS = 1;
	// To include multi-constellation support, grab the remaining values from rtklib.h
	static const unsigned int MAXSAT = NSATGPS;
	//static const unsigned int MAXSAT = (NSATGPS+NSATGLO+NSATGAL+NSATQZS+NSATCMP+NSATIRN+NSATSBS+NSATLEO);


	/** \brief RTKLIB earth parameters from ephemeris.c */
	static const double MU_GPS =	3.9860050e14;     	/* gravitational constant         ref [1] */
	static const double MU_GLO =	3.9860044e14;     	/* gravitational constant         ref [2] */
	static const double MU_GAL = 	3.986004418e14;   	/* earth gravitational constant   ref [7] */
	static const double MU_CMP =	3.986004418e14; 	/* earth gravitational constant   ref [9] */

	/** \brief RTKLIB earth parameters from ephemeris.c and rtklib.h */
	static const double OMGE_GPS =  7.2921151467e-5;	/* earth angular velocity (rad/s) */
	static const double OMGE_GLO = 	7.292115e-5;      	/* earth angular velocity (rad/s) ref [2] */
	static const double OMGE_GAL = 	7.2921151467e-5;  	/* earth angular velocity (rad/s) ref [7] */
	static const double OMGE_CMP = 	7.292115e-5; 		/* earth angular velocity (rad/s) ref [9] */


	/**
	 * The following data types are taken from RTKLIB
	 * to provide compatability with their tools and avoid reinventing the wheel
	 *
	 * Links:
	 * http://www.rtklib.com/
	 * https://github.com/tomojitakasu/RTKLIB
	 *
	 */


	// Taken from rtklib.h
	typedef struct {        /* time struct */
		time_t time;        /* time (s) expressed by standard time_t */
		double sec;         /* fraction of second under 1 s */
	} gtime_t;


	// Taken from rtklib.h
	typedef struct {        /* almanac type */
	    int sat;            /* satellite number */
	    int svh;            /* sv health (0:ok) */
	    int svconf;         /* as and sv config */
	    int week;           /* GPS/QZS: gps week, GAL: galileo week */
	    gtime_t toa;        /* Toa */
	                        /* SV orbit parameters */
	    double A,e,i0,OMG0,omg,M0,OMGd;
	    double toas;        /* Toa (s) in week */
	    double f0,f1;       /* SV clock parameters (af0,af1) */
	} alm_t;


	// Taken from rtklib.h
	typedef struct {        /* GPS/QZS/GAL broadcast ephemeris type */
		int sat;            /* satellite number */
		int iode,iodc;      /* IODE,IODC */
		int sva;            /* SV accuracy (URA index) */
		int svh;            /* SV health (0:ok) */
		int week;           /* GPS/QZS: gps week, GAL: galileo week */
		int code;           /* GPS/QZS: code on L2 */
							/* GAL: data source defined as rinex 3.03 */
							/* BDS: data source (0:unknown,1:B1I,2:B1Q,3:B2I,4:B2Q,5:B3I,6:B3Q) */
		int flag;           /* GPS/QZS: L2 P data flag */
							/* BDS: nav type (0:unknown,1:IGSO/MEO,2:GEO) */
		gtime_t toe,toc,ttr; /* Toe,Toc,T_trans */
							/* SV orbit parameters */
		double A,sqrt_A,e,e_sqr,i0,OMG0,omg,M0,deln,OMGd,idot;
		double crc,crs,cuc,cus,cic,cis;
		double toes, tocs;        /* Toe (s), Toc (s) in week */
		double fit;         /* fit interval (h) */
		double f0,f1,f2;    /* SV clock parameters (af0,af1,af2) */
		double tgd[4];      /* group delay parameters */
							/* GPS/QZS:tgd[0]=TGD */
							/* GAL    :tgd[0]=BGD E5a/E1,tgd[1]=BGD E5b/E1 */
							/* CMP    :tgd[0]=BGD1,tgd[1]=BGD2 */
		double Adot,ndot;   /* Adot,ndot for CNAV */
	} eph_t;


	typedef struct {
		double a0;
		double a1;
		double a2;
		double a3;

		double b0;
		double b1;
		double b2;
		double b3;
	} ionGPS_t;


	typedef struct {
		double A0;
		double A1;
		double T;
		double W;
	} utcGPS_t;


	struct ephSet_t {
		eph_t eph[CONST_PRN_MAX];				/* GPS ephemeris, stored with array index PRN-1 */
		bool ephValid[CONST_PRN_MAX] = {0};	/* 1 if the ephems of the corresponding eph are valid */
		double ephToes = -1;			/* The time of ephemerides for this set (redundant, but included for accessibility) */
										// TOE (s) should be int, but is kept double for RTKLIB compatibility
	    double utc_gps[4];  					/* GPS delta-UTC parameters {A0,A1,T,W} */
	    double ion_gps[8];  					/* GPS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
	    int leaps;          					/* leap seconds (s) */

	    // Deep-copies the source struct provided into this struct
	    void copyInto(ephSet_t src) {
	    	ephToes 	= src.ephToes;
	    	leaps 		= src.leaps;

	    	for (int i = 0; i < CONST_PRN_MAX; i++) {
	    		if (i < CONST_PRN_MAX) 	{ eph[i] = src.eph[i]; ephValid[i] = src.ephValid[i]; }
	    		if (i < 4) 				{ utc_gps[i] = src.utc_gps[i]; }
	    		if (i < 8) 				{ ion_gps[i] = src.ion_gps[i]; }
	    	}

	    	return;
	    }

	};
	typedef struct ephSet_t ephSet_t;


	typedef struct {
	    double utc_gps[4];  					/* GPS delta-UTC parameters {A0,A1,T,W} */
	    double ion_gps[8];  					/* GPS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
	    int leaps;
	} rinexHeaderParams_t;


	// Adapted from rtklib.h
	// A Lite version of the nav_t type so updating and transferring on GPU is cheaper
	struct navLite_t{        								/* navigation data type */
	    //eph_t eph[PRN_MAX];         	/* GPS ephemeris, stored with array index PRN-1 */
	    //bool eph_valid[PRN_MAX] = {0};/* 1 if the ephems of the corresponding eph are valid */
	    //utcGPS_t utc_gps;	  			/* GPS delta-UTC parameters {A0,A1,T,W} */
	    //ionGPS_t ion_gps;  			/* GPS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */

		// Consider moving the following into ephSet_t
		// (this way because RINEX only has one set of values for these params each day)
	    double utc_gps[4];  					/* GPS delta-UTC parameters {A0,A1,T,W} */
	    double ion_gps[8];  					/* GPS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
	    int leaps;          					/* leap seconds (s) */
	    ephSet_t *ephSet;						/* a set of ephemerides with the same toe */

	};
	typedef struct navLite_t navLite_t;



}
}

#endif

#ifndef __INC_RINEXPARSE_H_
#define __INC_RINEXPARSE_H_

#include "ephhelper.h"
#include "converters.h"
#include <stdio.h>

namespace dsp { namespace utils {

	int ReadRinex(FILE *fp, const char *opt, std::vector<dsp::utils::ephSet_t> &nav);


	// Make anonymous namespace to hide functions used by ReadRinex
	namespace {
		int ReadRinexHeader(FILE *fp, double *ver, char *type, int *sys, int *tsys, std::vector<dsp::utils::ephSet_t> &nav, dsp::utils::rinexHeaderParams_t *curParams);
		int DecodeRinexHeader(char *buff, dsp::utils::rinexHeaderParams_t *curParams);
		int ReadRinexBody(FILE *fp, const char *opt, double ver, int sys, std::vector<dsp::utils::ephSet_t> &nav, dsp::utils::rinexHeaderParams_t *curParams);
		int ParseRinexBody(FILE *fp, const char *opt, double ver, int sys, int *type, dsp::utils::eph_t *eph);
		int DecodeRinexBody(double ver, int sat, dsp::utils::gtime_t toc, const double *data, dsp::utils::eph_t *eph);

		int satno(int sys, int prn);
		int satsys(int sat, int *prn);
		static int uraindex(double value);
		int set_sysmask(const char *opt);
	}


	/** \brief RTKLIB RINEX length parameters */
	static const int MAXOBSTYPE = 64;					// Max number of obs type in RINEX
	static const int MAXRNXLEN 	= (16*MAXOBSTYPE+4);	// Max RINEX record length
	static const int EPH_SET_INC_SIZE = 16;				// How many sets of ephemerides to allocate space for
	static const int MAXPOSHEAD = 1024;					// Max head line position



}
}

#endif

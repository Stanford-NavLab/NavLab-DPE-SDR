#ifndef __INC_UTILS_STATEHELPER_H_
#define __INC_UTILS_STATEHELPER_H_

#include "consthelper.h"

namespace dsp { namespace utils {

	/** \brief 8D Position-Velocity-Time state (can initialize as double via State<double> for ex. or any other type)
	 *  Note: this is un-qualified (like statePT_t is) since it's the default type for DPE-related things.
	 *        statePT_t is qualified since it's "unusual" in DPE -- an abbreviated form for maninfolds. */
	template <typename T>
	struct state_t {
		T x;
		T y;
		T z;
		T delta_t;
		T x_dot;
		T y_dot;
		T z_dot;
		T delta_t_dot;
	};


	template <typename T>
	struct weightState_t {
		T a;
		T b;
		T c;
		T d;
		T score;
	};

	template <typename T>
	struct LLA_t {
		T lat;
		T lon;
		T alt;
	};


	typedef struct {
		state_t<double> *satPosPtr_d[CONST_PRN_MAX];
		bool valid_d[CONST_PRN_MAX] = {false}; // Used to mark if positions have been resolved for this
		int startTime_d; // TODO: Make this handle GPS week rollover (gtime_t not nice due to conversion needed in DPMeas)
		bool ready_d;
		bool replace_d;
	} satState_t;

}
}

#endif

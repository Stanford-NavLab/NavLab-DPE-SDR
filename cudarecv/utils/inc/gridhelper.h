#ifndef __INC_UTILS_GRIDHELPER_H_
#define __INC_UTILS_GRIDHELPER_H_

#include "consthelper.h"

namespace dsp { namespace utils {


	/** \brief 4D Position-Time state */
	template <typename T>
	struct statePosManifold_t {
		T x;
		T y;
		T z;
		T delta_t;
	};

	/** \brief 4D Velocity-TimeDrift state */
	template <typename T>
	struct stateVelManifold_t {
		T x_dot;
		T y_dot;
		T z_dot;
		T delta_t_dot;
	};


	/** \brief Enums for selecting different manifold grid spacing types
	 *         BatchCorrManifold will generate based on the value selected */
	enum ManifoldGridTypes
	{
		Uniform,
		Exponential,
		ArthurBasis
	};

}
}

#endif

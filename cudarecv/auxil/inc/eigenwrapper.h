
#ifndef INC__EIGENWRAPPER_H_
#define INC__EIGENWRAPPER_H_

#define EIGEN_NO_CUDA
#define EIGEN_NO_HALF
#include <Eigen/Core>

namespace auxil {

	template <class T>
	T* MakeIMatrix(int height, int width);

	template <class T>
	T* MakeZMatrix(int height, int width);

	double* MakeFMatrix(int height, int width, double SampleLength);

}


#endif

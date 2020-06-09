#include "eigenwrapper.h"



/**
 * Create an Identity matrix using Eigen and copy it
 * to the location specified by the input pointer
 *
 */
template <class T>
T* auxil::MakeIMatrix(int height, int width) {
	static Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> IMatWrap(height, width);
	IMatWrap = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(height, width);
	return IMatWrap.data();
	//memcpy(dataPtr, temp.data(), sizeof(T)*height*width);
}

/**
 * Create a Zero matrix using Eigen and copy it
 * to the location specified by the input pointer
 *
 */
template <class T>
T* auxil::MakeZMatrix(int height, int width) {
	static Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> ZMatWrap(height, width);
	ZMatWrap = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(height, width);
	return ZMatWrap.data();
	//memcpy(dataPtr, temp.data(), sizeof(T)*height*width);
}

/**
 * Create the F matrix using Eigen and copy it
 * to the location specified by the input pointer
 *
 */
double* auxil::MakeFMatrix(int height, int width, double SampleLength) {
	if (height >= 8 && width >= 8) {
		static Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> FMatWrap(height, width);
		FMatWrap = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>::Identity(height, width);
		FMatWrap.block<4,4>(0,4) = Eigen::Matrix<double,4,4>::Identity() * SampleLength;
		return FMatWrap.data();
	}
	else { return NULL; }

	//memcpy(dataPtr, temp.data(), sizeof(double)*height*width);
}


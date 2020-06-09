
#ifndef INC__CUEKF_H_
#define INC__CUEKF_H_

//#define EIGEN_NO_CUDA

#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
//#include <Eigen/Core>
#include <memory>
#include "cublas_v2.h"

#include "dsp.h"
#include "module.h"
//#include "naveng.h"
//#include "geocoord.h"
#include "constant.h"
//#include "eigenwrapper.h"

// Maximum EKF dimensionality
// (feel free to change this, but be advised performance may suffer)
#define EKF_MAX_DIM 16




// Constant for updating the Q matrix (move to DPMove later)
#define Q_CLOCK_DRIFT ((2.5e-10)*(2.5e-10)*WGS84_c*WGS84_c)

namespace dsp {
    //class cuEKF;
    //using Mat8 = Eigen::Matrix<double,8,8>;
    //using Vec8 = Eigen::Matrix<double,8,1>;

	// Based on http://www.shodor.org/media/content/petascale/materials/UPModules/matrixMultiplication/moduleDocument.pdf
	// Matrices stored in col-major order
	// M(row, col) = *(M.elements + col * M.height + row)
	struct cuMatrix {
		int width;
		int height;
		double** elements;

		cuMatrix() {
			width = 0;
			height = 0;
			elements = NULL;
		}

		cuMatrix(int initHeight, int initWidth, double** initElem) {
			width = initWidth;
			height = initHeight;
			elements = initElem;
		}
	};
	typedef struct cuMatrix cuMatrix;


};



class dsp::cuEKF : public dsp::Module {

	public:
		cuEKF(void);
		~cuEKF(void);

		int Start(void* cuFlowStream);

		int Stop(void);

		int Update(void* cuStream);

	protected:

		/**
		 * Functions
		 */

		// Run the steps of the KF
		int StepUpdate();
		int StepPredict();

		// These functions have the equations for f and h defined --
		// edit here to change motion or measurement model
		void GetfVal();
		void GethVal();

		// These functions get the values for u and Q at measurement k --
		// should eventually be hooked up to a module to process control inputs
		void GetuVal();
		int GetQVal();



		/**
		 * Variables
		 */

        /** Flag to signify if the call to Update() is the first call after Start() */
        char FirstUpdate = 1;
        cudaStream_t *cuStream;
        bool EnableEKF;


		// Default KF params for DPE
		//Eigen::Matrix<double,8,8> FDefDP;
		//Eigen::Matrix<double,8,8> HDefDP;
		//Eigen::Matrix<double,8,1> uDefDP;
		//Eigen::Matrix<double,8,8> QDefDP;
		//Eigen::Matrix<double,8,8> PDefDP;

		double SampleLength;

		dsp::fFunc_t fFun; // Pointer to function pointer for motion model

		char Started = 0;

		// GPU matrix operations
		cublasHandle_t handle;
		cublasStatus_t opStatus;

		// The current measurement
		int measIdx_k;
		// Store the previous measurement for "catch-up"
		int prev_measIdx_k;

		// Dimensions of this EKF
		int stateDim; // referred to in comments as "sd"
		int measDim;  // referred to in comments as "md"

		// State variables for this update
		dsp::cuMatrix X_k1k1;	// X_k-1|k-1 (sdx1)
		double* X_k1k1_Arr_d;
		dsp::cuMatrix X_kk1; 	// X_k|k-1 (sdx1)
		//double **X_kk1_Arr_Ptr;
		double* X_kk1_Arr_d;
		dsp::cuMatrix X_kk;		// X_k|k (sdx1)
		double* X_kk_Arr_d;

		// Covariance variables for this update
		dsp::cuMatrix P_k1k1;	// P_k-1|k-1 (sdxsd)
		double* P_k1k1_Arr_d;
		dsp::cuMatrix P_kk1;	// P_k|k-1 (sdxsd)
		//double **P_kk1_Arr_Ptr;
		double* P_kk1_Arr_d;
		dsp::cuMatrix P_kk;		// P_k|k (sdxsd)
		double* P_kk_Arr_d;
		dsp::cuMatrix P_temp;	// Temporary variable to hold I
		double* P_temp_Arr_d;


		// Variables for DPMeas
		dsp::cuMatrix z; // (mdx1)
		double *zPtr_d;
		//double* z_Arr_d; 	// Don't need z_arr since it comes from the pipeline
		dsp::cuMatrix R; // (mdxmd)
		double *RPtr_d;
		//double* R_Arr_d;	// Don't need R_Arr since it comes from the pipeline

		// Variables for DPMove
		dsp::cuMatrix u; // (sdx1)
		double* u_Arr_d;
		dsp::cuMatrix Q; // (sdxsd)
		double* Q_Arr_d;
		dsp::cuMatrix Q_temp; // (sdxsd)
		double* Q_temp_Arr_d;

		// Variables for DPModel
		dsp::cuMatrix fval; // (sdx1)
		double* fval_Arr_d;
		dsp::cuMatrix F; // (sdxsd)
		double* F_Arr_d;
		dsp::cuMatrix F_temp;	// Temporary variable to hold F_k*P_k-1|k-1
		double* F_temp_Arr_d;
		dsp::cuMatrix hval; // (mdx1)
		double* hval_Arr_d;
		dsp::cuMatrix H; // (mdxmd)
		double* H_Arr_d;


		// Update step variables
		dsp::cuMatrix y; // (mdx1)
		double* y_Arr_d;
		dsp::cuMatrix S; // (mdxmd)
		double* S_Arr_d;
		double* S_DevArr_d[1]; // Device array of device pointer to S (for ri/rf in inverse)
		double** S_DevArrPtr_d; // Host pointer to device array of device pointer to S (for ri/rf in inverse)
		dsp::cuMatrix S_inv; // (mdxmd)
		double* S_inv_Arr_d;
		double* S_inv_DevArr_d[1]; // Device array of device pointer to S_inv (for ri/rf in inverse)
		double** S_inv_DevArrPtr_d; // Host pointer to device array of device pointer to S_inv (for ri/rf in inverse)
		dsp::cuMatrix S_temp; // Temporary variable to hold H_k*P_k|k-1
		double* S_temp_Arr_d;
		int* S_rf_info_d;
		int S_rf_info_host[1];
		int* rf_pivot_d;
		int* S_ri_info_d;
		int S_ri_info_host[1];
		int* ri_pivot_d;
		dsp::cuMatrix K; // (sdxmd)
		double* K_Arr_d;
		dsp::cuMatrix K_temp;
		double* K_temp_Arr_d;


		// Other matrices
		dsp::cuMatrix I_P; // Identity matrix of the same size as P (sdxsd)
		double* I_P_Arr_d;


		// Store the array sizes
		int size_sdx1;
		int size_sdxsd;
		int size_mdx1;
		int size_mdxmd;



		// Because cuBLAS is silly, it wants pointers to scalar values used in gemm and the like
		const double constNegOne = -1.0;
		const double constZero = 0.0;
		const double constPosOne = 1.0;



		// Move to DPMove later -- for updating Q
		double *upQ_LPFv_vals_d;
		int *upQ_LPFv_idx_d;
		double *upQ_LPFv_avg_d;


};

#endif

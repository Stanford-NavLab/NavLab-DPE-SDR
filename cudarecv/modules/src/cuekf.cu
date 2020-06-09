
//#include <cuda_runtime.h>
//#include <cuda_runtime_api.h>
#include "cuekf.h"
#include <cstring>
#include <iostream>
//#include "cublas_v2.h"
#include "errorhandler.h"
#include "eigenwrapper.cpp" // Keep this included as .cpp, not as .h!



/**
 * ==============================================================
 * || MODULE SUMMARY                                           ||
 * ==============================================================
 *
 * Convention:
 *    k1k1: time k-1|k-1 corresponds to the outcome of the previous measurement step (input to this step)
 *    kk1:	time k|k-1	 corresponds to the outcome of the current predict step
 *    kk:	time k|k	 corresponds to the outcome of the current measurement step
 *
 *
 * Idea: Create three types of modules -- DPMeas, DPModel, and DPMove
 * 	  DPMeas generates the measurement -- z and R (implemented as BCS and BCM)
 * 	  DPModel model's the receiver's motion and measurements -- f and h (not implemented)
 * 	  DPMove captures the receiver's inputs to the world -- u and Q (not implemented)
 * Hooking these up to the Kalman filter gives a robot action flow.
 * Some work has begun in this direction, but the authors have thus far focused on the DPE algorithm -- DPMeas only
 *
 */



/**
 * \brief Device-side update to the Q matrix
 *
 *  This algorithm comes from www.mathworks.com/help/control/examples/state-estimation-using-time-varying-kalman-filter.html
 *  and is implemented in PyGNSS since Yuting Ng created it for a multi-receiver vector tracking implementation
 *  (https://web.stanford.edu/~gracegao/publications//journal//2017%20IEEE%20TAES_MRVT.pdf)
 */
__global__ void
EKF_Update_Q(double *LPF_vals, int *LPF_idx, double *LPF_avg, double *Q_vals, int Q_dim, double *X_ECEF) {

	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int stride = blockDim.x * gridDim.x;

    __shared__ double LPF_v_val;

    // Compute the value of v
    if (i == 0) {
    	// Compute the average
    	double v = norm(3, X_ECEF+4);
    	*LPF_avg = *LPF_avg - (LPF_vals[*LPF_idx]) + (v/20.0);
    	LPF_vals[*LPF_idx] = v/20.0;
    	++(*LPF_idx); // Increment the LPF index value (by first dereferencing pointer)
    	if ((*LPF_idx) >= 20) { (*LPF_idx) = 0; }
    	// Compute the value of v
    	LPF_v_val = 1.0 + 250.0/(fmin(fmax((*LPF_avg) * (*LPF_avg), 50.0),125.0));

    }

    // Hold here until thread 0 has calculated the new value of v
    __syncthreads();

    // Change the corresponding values of Q
    while (i < (Q_dim * Q_dim)) {
    	if (i == (4*8 + 4) || i == (5*8 + 5) || i == (6*8 + 6)) {
    		Q_vals[i] = LPF_v_val;
    	}
    	else if (i == 63) {
    		Q_vals[i] = Q_CLOCK_DRIFT;
    	}
    	else {
    		Q_vals[i] = 0;
    	}

    	i += stride;
    }

}


/**
 * \brief Device-side identity matrix generator
 */
__global__ void
EKF_MakeIMatrix(double *arr, int arrHeight) {

	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int stride = blockDim.x * gridDim.x;

    while (i < arrHeight*arrHeight) {

    	arr[i] = 0;
    	for (int j = 0; j < arrHeight; j++) {
    		if (i == j*(arrHeight+1)) {
    			arr[i] = 1;
    			continue;
    		}
    	}

    	i += stride;
    }
}


/**
 * \brief Device-side generator a random walk matrix for a state with 4 position-time dimensions and 4 velocity-drift dimensions.
 */
__global__ void
EKF_MakeDPERandomWalkFMatrix(double *arr, double samplePeriod) {

	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int stride = blockDim.x * gridDim.x;

    // This is a specific matrix, so size values are hardcoded
    while (i < 64) {

    	arr[i] = 0;
    	for (int j = 0; j < 8; j++) {
    		if (i == j*(8+1)) {
    			arr[i] = 1;
    			continue;
    		}
    	}

    	for (int j = 0; j < 4; j++) {
    		// T corresponds to the dX-d_deltaT elements of the X-delta_t states,
    		// so, since cuBLAS is column-major, place them accordingly as below
    		if (i == j*(8+1)+32) {
    			arr[i] = samplePeriod;
    			continue;
    		}
    	}

    	i += stride;
    }
}


/**
 * \brief Device-side function to pass the measurement from the last iteration to the next iteration
 *
 * (When you want to disable the EKF)
 */
__global__ void
EKF_PassMeas(double *meas, double *x_k1k1, double *x_kk1, int numStates) {

	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int stride = blockDim.x * gridDim.x;

    while (i < numStates*2) {
    	x_k1k1[i] = meas[i];
    	x_kk1[i] = meas[i];

    	i += stride;
    }
}


/** \brief Debug function to make viewing device memory easier. Placing a breakpoint in a kernel doesn't guarantee thread safety.
 * 		   This function makes no modifications to data, so run it on a stream to view device memory at that point.
 *
 * \param *arr		The double array you want to view
 * \param arrLen	The size of the array you want to view
 * \param anc		Give this a unique value to identify where you are in the stream
 *
 */
__global__ void
EKF_valInspector(double *arr, int arrLen, int anc)
{
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int stride = blockDim.x * gridDim.x;

    double test2 = 0;

    while (i<arrLen) {

		double test = arr[i];

		if (i == 0) {
			test2 = 1;
		}

		i += stride;
    }

}

/** \brief Debug function to make viewing device memory easier. Placing a breakpoint in a kernel doesn't guarantee thread safety.
 * 		   This function makes no modifications to data, so run it on a stream to view device memory at that point.
 *
 * \param **arr		The double-pointed array you want to view
 * \param arrLen	The size of the array you want to view
 * \param anc		Give this a unique value to identify where you are in the stream
 *
 */
__global__ void
EKF_valInspector2(double **arr, int arrLen, int anc)
{
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int stride = blockDim.x * gridDim.x;

    //int currChan;
    //int currSamp;

    double test2 = 0;

    while (i<arrLen) {

        double* test = arr[i];

        if (i == 0) {
        	test2 = 1;
        }

    	i += stride;
    }

}



dsp::cuEKF::cuEKF() {

    ModuleName = "cuEKF";
    AllocateInputs(9);
    AllocateOutputs(3);

    Started = 0;

    /**
     * INPUT CONFIGURATION
     */

    // From DPInit: the initial state x and covariance P at a designed measurement index k
    ConfigExpectedInput(0, "InitX", DOUBLE_t, STATE, VECTORLENGTH_ANY);
    ConfigExpectedInput(1, "InitP", DOUBLE_t, COVARIANCE, VECTORLENGTH_ANY);
    ConfigExpectedInput(2, "InitK", INT_t, VALUE, 1);

    // From DPMeas: the current measurement z_k and its covariance R_k
    ConfigExpectedInput(3, "zVal", DOUBLE_t, STATE, VECTORLENGTH_ANY);
    ConfigExpectedInput(4, "RVal", DOUBLE_t, COVARIANCE, VECTORLENGTH_ANY); // See how PyGNSS computes R!!!!!!!

    // From DPMove: pointers to functions to return the current movement u and its covariance Q
    // Functions take:
    //		measurement counter int k
    //		matrix u location double* elements (to update)
    //  	matrix size int dim
    //ConfigExpectedInput(5, "uFunc", DOUBLE_t, FUNCTION_PTR, 1);
    //ConfigExpectedInput(6, "QFunc", DOUBLE_t, FUNCTION_PTR, 1);

    // From DPModel: pointers to functions to return values of functions movement model f
    // and measurement model h
    // fFunc takes:
    // 		matrix x_k-1|k-1 location double* elements
    //		matrix u_k location double* elements
    //		matrix size int dim
    //		matrix fval location double* elements (to update)
    // 		matrix F location double* elements (to update)
    //ConfigExpectedInput(7, "fFunc", DOUBLE_t, FUNCTION_PTR, 1);
    // hFunc takes:
    //		matrix x_k|k-1 location double* elements
    //		matrix size int dim
    //		matrix hval location double* elements (to update)
    // 		matrix H location double* elements (to update)
    //ConfigExpectedInput(8, "hFunc", DOUBLE_t, FUNCTION_PTR, 1);



    /**
     * OUTPUT CONFIGURATION
     */

    // The current state x_k|k-1
    ConfigOutput(0, "xCurrkk1", DOUBLE_t, STATE, CUDA_DEVICE, VECTORLENGTH_ANY, NULL, 0);
    ConfigOutput(1, "PCurrkk1", DOUBLE_t, COVARIANCE, CUDA_DEVICE, VECTORLENGTH_ANY, NULL, 0);
    ConfigOutput(2, "xCurrk1k1", DOUBLE_t, STATE, CUDA_DEVICE, VECTORLENGTH_ANY, NULL, 0);


    /**
     *  PARAMETER CONFIGURATION
     */
    // The samplelength
    InsertParam("SampleLength", (void*)&SampleLength, DOUBLE_t, sizeof(double), sizeof(double));
    InsertParam("EnableEKF",	(void*)&EnableEKF,		BOOL_t, sizeof(bool), sizeof(bool));

}

dsp::cuEKF::~cuEKF() {
	if (Started) Stop();
    delete [] inputs;
    delete [] outputs;
    delete [] expectedInputs;
}


int dsp::cuEKF::Start(void* cuFlowStream) {
	// Check module status and report accordingly
	if (Started) {
        std::clog << "[" << ModuleName << "] Start: Already Started." << std::endl;
        return 0;
    }
    std::clog << "[" << ModuleName << "] Starting ... " << std::flush;


	// Set the CUDA Stream for the GPU operations
    cuStream = (cudaStream_t*)cuFlowStream;


    // Figure out how big the filter will be
    stateDim 	= inputs[0]->VectorLength;
    measDim		= inputs[3]->VectorLength;

    // Don't let the user initialize a filter size that's too large for computation purposes on the NVIDIA Jetson TX2
    if(stateDim > EKF_MAX_DIM || measDim > EKF_MAX_DIM) {
    	std::clog << "[" << ModuleName << "] Your desired filter size has " << stateDim <<
    			" and " << measDim <<
    			" dimensionality! To ensure computation doesn't blow up, the filter size" <<
    			" is capped at dimensionality " << EKF_MAX_DIM << ". If you're sure you want" <<
    			" your size, change EKF_MAX_DIM in modules/cuekf.h." << std::endl;
    	return -1;
    }




    // Create the cublas handle for matrix operations
    cuBLASCheckMSt(cublasCreate(&handle));


    // Initialize CUDA device memory for the filter

    // Initialize sdx1 vectors
    size_t size = sizeof(double)*stateDim;
    size_sdx1 = size;
    cuCheckMSt(cudaMalloc((void**)&X_k1k1_Arr_d, size));
    cuCheckMSt(cudaMemcpyAsync(X_k1k1_Arr_d, (double*)inputs[0]->Data, size, cudaMemcpyHostToDevice, *cuStream));
    //cuCheckMSt(cudaMallocManaged((void**)&X_kk1_Arr_d, size)); // Consider for GPU->CPU transfers in logging
    cuCheckMSt(cudaMalloc((void**)&X_kk1_Arr_d, size));
    cuCheckMSt(cudaMemcpyAsync(X_kk1_Arr_d, (double*)inputs[0]->Data, size, cudaMemcpyHostToDevice, *cuStream));
    cuCheckMSt(cudaMalloc((void**)&X_kk_Arr_d, size));
    cuCheckMSt(cudaMemcpyAsync(X_kk_Arr_d, (double*)inputs[0]->Data, size, cudaMemcpyHostToDevice, *cuStream));
    cuCheckMSt(cudaMalloc((void**)&u_Arr_d, size));
    cuCheckMSt(cudaMalloc((void**)&fval_Arr_d, size));

    // Initialize sdxsd matrices
    size = size * stateDim;
    size_sdxsd = size;
    cuCheckMSt(cudaMalloc((void**)&P_k1k1_Arr_d, size));
    cuCheckMSt(cudaMemcpy(P_k1k1_Arr_d, (double*)inputs[1]->Data, size, cudaMemcpyHostToDevice));
    cuCheckMSt(cudaMalloc((void**)&P_kk1_Arr_d, size));
    cuCheckMSt(cudaMalloc((void**)&P_kk_Arr_d, size));
    cuCheckMSt(cudaMalloc((void**)&P_temp_Arr_d, size));
    cuCheckMSt(cudaMalloc((void**)&Q_Arr_d, size));
    cuCheckMSt(cudaMalloc((void**)&Q_temp_Arr_d, size));
    cuCheckMSt(cudaMalloc((void**)&F_Arr_d, size));
    cuCheckMSt(cudaMalloc((void**)&F_temp_Arr_d, size));
    cuCheckMSt(cudaMalloc((void**)&I_P_Arr_d, size));

    // Initialize mdx1 vectors
    size = sizeof(double)*measDim;
    size_mdx1 = size;
    //cuCheckMSt(cudaMalloc((void**)&z_Arr_d, size)); // Don't need z since it comes from the pipeline
    cuCheckMSt(cudaMalloc((void**)&hval_Arr_d, size));
    cuCheckMSt(cudaMalloc((void**)&y_Arr_d, size));
    cuCheckMSt(cudaMalloc((void**)&rf_pivot_d, sizeof(int)*measDim));
    cuCheckMSt(cudaMalloc((void**)&ri_pivot_d, sizeof(int)*measDim));

    // Initialize mdxmd vectors
    size = size * measDim;
    size_mdxmd = size;
    //cuCheckMSt(cudaMalloc((void**)&R_Arr_d, size)); // Don't need R since it comes from the pipeline
    cuCheckMSt(cudaMalloc((void**)&H_Arr_d, size));
	cuCheckMSt(cudaMalloc((void**)&S_Arr_d, size));
	cuCheckMSt(cudaMalloc((void**)&S_inv_Arr_d, size));
	cuCheckMSt(cudaMalloc((void**)&S_temp_Arr_d, size));
	// Store the pointer to S/S_inv using a host pointer to a device array to a device pointer to S/S_inv
	// (per ri/rf inverse specifications, see https://devblogs.nvidia.com/cuda-pro-tip-how-call-batched-cublas-routines-cuda-fortran/)
	// though this implementation is only one EKF, so there's only one element allocated in the array
	S_DevArr_d[0] = S_Arr_d;
	cuCheckMSt(cudaMalloc((void**)&S_DevArrPtr_d, sizeof(S_DevArr_d)));
	cuCheckMSt(cudaMemcpy(S_DevArrPtr_d, S_DevArr_d, sizeof(S_DevArr_d), cudaMemcpyHostToDevice));
	S_inv_DevArr_d[0] = S_inv_Arr_d;
	cuCheckMSt(cudaMalloc((void**)&S_inv_DevArrPtr_d, sizeof(S_inv_DevArr_d)));
	cuCheckMSt(cudaMemcpy(S_inv_DevArrPtr_d, S_inv_DevArr_d, sizeof(S_inv_DevArr_d), cudaMemcpyHostToDevice));

	// Initialize sdxmd vectors
	size = sizeof(double)*stateDim*measDim;
    cuCheckMSt(cudaMalloc((void**)&K_Arr_d, size));
    cuCheckMSt(cudaMalloc((void**)&K_temp_Arr_d, size));


    // Debug
    cuCheckMSt(cudaMalloc((void**)&S_rf_info_d, sizeof(int)));
    cuCheckMSt(cudaMalloc((void**)&S_ri_info_d, sizeof(int)));


    // Make sure all allocation operations have completed before making the matrices
    cuCheckMSt(cudaStreamSynchronize(*cuStream));
    cuCheckMSt(cudaDeviceSynchronize());


    // Create matrices for all these variables

    // Create the sdx1 vectors
    X_k1k1 = dsp::cuMatrix(stateDim, 1, &X_k1k1_Arr_d);
    X_kk1 = dsp::cuMatrix(stateDim, 1, &X_kk1_Arr_d);
    X_kk = dsp::cuMatrix(stateDim, 1, &X_kk_Arr_d);
    u = dsp::cuMatrix(stateDim, 1, &u_Arr_d);
    fval = dsp::cuMatrix(stateDim, 1, &fval_Arr_d);

    // Create the sdxsd matrices
    P_k1k1 = dsp::cuMatrix(stateDim, stateDim, &P_k1k1_Arr_d);
    P_kk1 = dsp::cuMatrix(stateDim, stateDim, &P_kk1_Arr_d);
    P_kk = dsp::cuMatrix(stateDim, stateDim, &P_kk_Arr_d);
    P_temp = dsp::cuMatrix(stateDim, stateDim, &P_temp_Arr_d);
    Q = dsp::cuMatrix(stateDim, stateDim, &Q_Arr_d);
    Q_temp = dsp::cuMatrix(stateDim, stateDim, &Q_temp_Arr_d);
    F = dsp::cuMatrix(stateDim, stateDim, &F_Arr_d);
    F_temp = dsp::cuMatrix(stateDim, stateDim, &F_temp_Arr_d);
    I_P = dsp::cuMatrix(stateDim, stateDim, &I_P_Arr_d);

    // Create the mdx1 vectors
    zPtr_d = (double*)inputs[3]->Data;
    z = dsp::cuMatrix(measDim, 1, &zPtr_d);
    hval = dsp::cuMatrix(measDim, 1, &hval_Arr_d);
    y = dsp::cuMatrix(measDim, 1, &y_Arr_d);

    // Create the mdxmd matrices
    RPtr_d = (double*)inputs[4]->Data;
    R = dsp::cuMatrix(measDim, measDim, &RPtr_d);
    H = dsp::cuMatrix(measDim, measDim, &H_Arr_d);
    S = dsp::cuMatrix(measDim, measDim, &S_Arr_d);
    S_inv = dsp::cuMatrix(measDim, measDim, &S_inv_Arr_d);
    S_temp = dsp::cuMatrix(measDim, measDim, &S_temp_Arr_d);

    // Create the sdxmd matrices
    K = dsp::cuMatrix(stateDim, measDim, &K_Arr_d);
    K_temp = dsp::cuMatrix(stateDim, measDim, &K_temp_Arr_d);



    // Get the other variables
    measIdx_k = 0; // Default value until this can be passed in -- other modules must be initialized first
    prev_measIdx_k = -1;

    // The cool way to do motion model updates is to get a pointer to function pointer (not implemented yet)
    //fFun = (dsp::fFunc_t*)inputs[3]->Data;


    // Initialize other matrices
    cuCheckMSt(cudaMemcpyAsync(I_P_Arr_d, auxil::MakeIMatrix<double>(I_P.height, I_P.width),
    		sizeof(double)*stateDim*stateDim, cudaMemcpyHostToDevice, *cuStream));
    EKF_MakeIMatrix<<<1, 64, 0, *cuStream>>>(*I_P.elements, I_P.height);


    // Create default DPE KF matrices
    EKF_MakeDPERandomWalkFMatrix<<<1, 64, 0, *cuStream>>>(*F.elements, SampleLength);
    EKF_MakeIMatrix<<<1, 64, 0, *cuStream>>>(*H.elements, H.height);
    cuCheckMSt(cudaMemset((void*)(*u.elements), 0, 	sizeof(double)*u.height*u.width));
    EKF_MakeIMatrix<<<1, 64, 0, *cuStream>>>(*Q.elements, Q.height);
    EKF_MakeIMatrix<<<1, 64, 0, *cuStream>>>(*P_kk1.elements, P_kk1.height);


    // Wait for these matrices to be built
    cuCheckMSt(cudaDeviceSynchronize());


    // Set up Q update (this should be moved to DPMove when DPMove is implemented)
    cuCheckMSt(cudaMalloc((void**)&upQ_LPFv_vals_d, 	sizeof(double)*20));
    cuCheckMSt(cudaMalloc((void**)&upQ_LPFv_idx_d, 		sizeof(int)));
    cuCheckMSt(cudaMalloc((void**)&upQ_LPFv_avg_d, 		sizeof(double)));
    cuCheckMSt(cudaMemset((void*)upQ_LPFv_vals_d, 	0, 	sizeof(double)*20));
    cuCheckMSt(cudaMemset((void*)upQ_LPFv_idx_d, 	0, 	sizeof(int)));
    cuCheckMSt(cudaMemset((void*)upQ_LPFv_avg_d, 	0, 	sizeof(double)));


    // Update outputs now that space has been allocated
    outputs[0].Data = (void*)X_kk1_Arr_d;
    outputs[0].VectorLength = X_kk1.width * X_kk1.height;
    outputs[1].Data = *P_kk1.elements;
    outputs[1].VectorLength = P_kk1.width * P_kk1.height;
    outputs[2].Data = *X_k1k1.elements;
    outputs[2].VectorLength = X_k1k1.width * X_k1k1.height;

    // Make sure all GPU tasks have completed before continuing
    cuCheckMSt(cudaDeviceSynchronize());


    // The initial value is a k-1|k-1 measurement; predict the next state before continuing
    if (EnableEKF) {
		int stepPredictRetVal = StepPredict();
		if (stepPredictRetVal != 0) { return -1; }
    }


    // Make sure all GPU tasks have completed before continuing
    cuCheckMSt(cudaStreamSynchronize(*cuStream));
    cuCheckMSt(cudaDeviceSynchronize());

    // Signifies that the next call to update() will be the first after start()
    Started = 1;

    std::clog << "Started." << std::endl;
    return 0;

}

int dsp::cuEKF::Stop() {
    int ret = 0;
    if (Started == 0) {
        std::clog << "[" << ModuleName << "] Stop: Wasn't Started." << std::endl;
        return 0;
    }
    std::clog << "[" << ModuleName << "] Stopping ... " << std::flush;

    // Free device memory
    cuCheckMSp(cudaFree((void*)X_k1k1_Arr_d));
    cuCheckMSp(cudaFree((void*)X_kk1_Arr_d));
    //cuCheckMSp(cudaFree((void*)X_kk_Arr_d));
    cuCheckMSp(cudaFree((void*)P_k1k1_Arr_d));
    cuCheckMSp(cudaFree((void*)P_kk1_Arr_d));
    //cuCheckMSp(cudaFree((void*)P_kk_Arr_d));
    cuCheckMSp(cudaFree((void*)P_temp_Arr_d));
    cuCheckMSp(cudaFree((void*)u_Arr_d));
    cuCheckMSp(cudaFree((void*)Q_Arr_d));
    cuCheckMSp(cudaFree((void*)Q_temp_Arr_d));
    cuCheckMSp(cudaFree((void*)fval_Arr_d));
    cuCheckMSp(cudaFree((void*)F_Arr_d));
    cuCheckMSp(cudaFree((void*)F_temp_Arr_d));
    cuCheckMSp(cudaFree((void*)hval_Arr_d));
    cuCheckMSp(cudaFree((void*)H_Arr_d));
    cuCheckMSp(cudaFree((void*)y_Arr_d));
    cuCheckMSp(cudaFree((void*)S_Arr_d));
    cuCheckMSp(cudaFree((void*)S_temp_Arr_d));
    cuCheckMSp(cudaFree((void*)S_inv_Arr_d));
    cuCheckMSp(cudaFree((void*)K_Arr_d));
    cuCheckMSp(cudaFree((void*)K_temp_Arr_d));
    cuCheckMSp(cudaFree((void*)I_P_Arr_d));
    cuCheckMSp(cudaFree((void*)upQ_LPFv_vals_d));
    cuCheckMSp(cudaFree((void*)upQ_LPFv_idx_d));
    cuCheckMSp(cudaFree((void*)upQ_LPFv_avg_d));
    cuCheckMSp(cudaFree((void*)rf_pivot_d));
    cuCheckMSp(cudaFree((void*)ri_pivot_d));
    cuCheckMSp(cudaFree((void*)S_DevArrPtr_d));
    cuCheckMSp(cudaFree((void*)S_inv_DevArrPtr_d));


    cublasDestroy(handle);

    Started = 0;
    std::clog << "Stopped." << std::endl;

    return ret;
}

int
dsp::cuEKF::Update(void* cuFlowStream) {

    if (Started == 0){
        std::cerr << "[" << ModuleName
                  << "] Error: Update() Failed due to EKF not initialized"
                  << std::endl;
        return -1;
    }

    if (FirstUpdate) {
        cuStream = (cudaStream_t*)cuFlowStream;
    	cuBLASCheckMSt(cublasSetStream(handle, *cuStream)); // put this on a specific stream??
    	FirstUpdate = 0;
    }


    // TODO: put kernel calls on the cuStream since all modules further down require it???
    if (EnableEKF) {
		// Run KF
		int stepUpdateRetVal = StepUpdate();
		if (stepUpdateRetVal != 0) { return -1;	}

		while (prev_measIdx_k < measIdx_k) {
			int stepPredictRetVal = StepPredict();
			if (stepPredictRetVal != 0) { return -1; }
			// Increment the measurement counter
			prev_measIdx_k += 1;
		}
		measIdx_k += 1;
    }
    else {
    	EKF_PassMeas<<<1, 32, 0, *cuStream>>>(*z.elements, *X_k1k1.elements, *X_kk1.elements, 8);
    }


    cuCheckMSt(cudaStreamSynchronize(*cuStream));


    return 0;
}

/**
 * This function would contain the mathematical equations for the motion model
 * (state transition equations)
 *
 * This would compute x_k|k-1 given x_k-1|k-1 and u_k and store it in fval
 * as well as the Jacobian at k and store it in F
 *
 * To upgrade, for accessibility, the transition algorithms should be
 * defined in another module and accessed via a function pointer
 * or similar methodology.
 *
 */
void
dsp::cuEKF::GetfVal() {

	// Assume random walk state transition model (get this from a function later)



	return;
}


// Computes k|k-1 from k-1|k-1 (last step's k|k)
int
dsp::cuEKF::StepPredict() {

	// Update Q
	GetQVal();
    cuCheckMSt(cudaStreamSynchronize(*cuStream));

	// Cool way to do updates: call function passed as a pointer to function pointer
	//(*fFun)(X_k1k1.elements, u.elements, fval.elements, F.elements);

	// Compute x_k|k-1
	cuBLASCheckMSt(cublasDgemv(handle, CUBLAS_OP_N, stateDim, stateDim,
			&constPosOne, *F.elements, F.height, *X_k1k1.elements, 1, &constZero, *X_kk1.elements, 1));

	// Compute F*P_k-1|k-1 then multiply that by F^T + Q_k
	cuCheckMSt(cudaMemcpy(*P_kk1.elements, *Q.elements, size_sdxsd, cudaMemcpyDeviceToDevice));
	cuBLASCheckMSt(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, F.height, P_k1k1.width, F.height,
				&constPosOne, *F.elements, F.height, *P_k1k1.elements, P_k1k1.height, &constZero, *F_temp.elements, F_temp.height));
	cuBLASCheckMSt(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_T, F_temp.height, F.width, F_temp.height,
			&constPosOne, *F_temp.elements, F_temp.height, *F.elements, F.height, &constPosOne, *P_kk1.elements, P_kk1.height));

    cuCheckMSt(cudaStreamSynchronize(*cuStream));


	// Debug: Make Q identity
    //EKF_MakeIMatrix<<<1, 64, 0, *cuStream>>>(*Q.elements, Q.height);

	//EKF_valInspector<<<1, 128, 0, *cuStream>>>(*Q.elements, 8, 43);


	return 0;
}

// Computes k|k from k|k-1
int
dsp::cuEKF::StepUpdate() {

	// Compute y_k
	cuCheckMSt(cudaMemcpy(*y.elements, *z.elements, size_mdx1, cudaMemcpyDeviceToDevice));
	cuBLASCheckMSt(cublasDgemv(handle, CUBLAS_OP_N, H.height, H.width,
			&constNegOne, *H.elements, H.height, *X_kk1.elements, 1, &constPosOne, *y.elements, 1));


	// Compute S_k
	cuCheckMSt(cudaMemcpy(*S.elements, *R.elements, size_mdxmd, cudaMemcpyDeviceToDevice));
	cuBLASCheckMSt(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, H.height, P_kk1.width, H.width,
			&constPosOne, *H.elements, H.height, *P_kk1.elements, P_kk1.height, &constZero,* S_temp.elements, S_temp.height));
	cuBLASCheckMSt(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_T, S_temp.height, H.width, S_temp.width,
			&constPosOne, *S_temp.elements, S_temp.height, *H.elements, H.height, &constPosOne, *S.elements, S.height));


	// Compute K
	// Compute S^-1
	// TODO: info is passed as a pointer; need to copy from device memory to host to be read
	cuBLASCheckMSt(cublasDgetrfBatched(handle, S.height, S_DevArrPtr_d, S.height, rf_pivot_d, S_rf_info_d, 1));
	cuCheckMSt(cudaMemcpy(S_rf_info_host,S_rf_info_d,sizeof(int),cudaMemcpyDeviceToHost));
	if (S_rf_info_host[0] != 0) {
        std::cerr << "[" << ModuleName
                  << "] Error: StepUpdate() S inversion failed at rf with error " << S_rf_info_host[0]
                  << std::endl;
        return -1;
	}
	cuBLASCheckMSt(cublasDgetriBatched(handle, S.height, (const double **)S_DevArrPtr_d, S.height, rf_pivot_d, S_inv_DevArrPtr_d, S_inv.height, S_ri_info_d, 1));
	cuCheckMSt(cudaMemcpy(S_ri_info_host,S_ri_info_d,sizeof(int),cudaMemcpyDeviceToHost));
	if (S_ri_info_host[0] != 0) {
        std::cerr << "[" << ModuleName
                  << "] Error: StepUpdate() S inversion failed at ri with error " << S_ri_info_host[0]
                  << std::endl;
        return -1;
	}


	// Compute P_k|k-1*H_k^T then multiply by S_inv
	cuBLASCheckMSt(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_T, P_kk1.height, H.width, P_kk1.width,
			&constPosOne, *P_kk1.elements, P_kk1.height, *H.elements, H.height, &constZero, *K_temp.elements, K_temp.height));
	cuBLASCheckMSt(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, K_temp.height, S_inv.width, K_temp.width,
			&constPosOne, *K_temp.elements, K_temp.height, *S_inv.elements, S_inv.height, &constZero, *K.elements, K.height));


	// Compute x_k|k = x_k-1|k-1 for the predict step
	cuCheckMSt(cudaMemcpy(*X_k1k1.elements, *X_kk1.elements, size_sdx1, cudaMemcpyDeviceToDevice));
	cuBLASCheckMSt(cublasDgemv(handle, CUBLAS_OP_N, K.height, K.width,
				&constPosOne, *K.elements, K.height, *y.elements, 1, &constPosOne, *X_k1k1.elements, 1));


	// Compute P_k|k = P_k-1|k-1 for the predict step
	cuCheckMSt(cudaMemcpy(*P_temp.elements, *I_P.elements, size_sdxsd, cudaMemcpyDeviceToDevice));
	cuBLASCheckMSt(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, K.height, H.width, K.width,
			&constNegOne, *K.elements, K.height, *H.elements, H.height, &constPosOne, *P_temp.elements, P_temp.height));
	cuBLASCheckMSt(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, P_temp.height, P_kk1.width, P_temp.width,
			&constPosOne, *P_temp.elements, P_temp.height, *P_kk1.elements, P_kk1.height, &constZero, *P_k1k1.elements, P_k1k1.height));


    cuCheckMSt(cudaStreamSynchronize(*cuStream));

	return 0;
}


/** The functions in here should be moved to DPMeas and replaced with a simple function call
 *  (or something similar in effect) to the DPMove module.
 *
 *  This algorithm comes from www.mathworks.com/help/control/examples/state-estimation-using-time-varying-kalman-filter.html
 *  and is implemented in PyGNSS since Yuting Ng created it for a multi-receiver vector tracking implementation
 *  (https://web.stanford.edu/~gracegao/publications//journal//2017%20IEEE%20TAES_MRVT.pdf)
 *
 */
int
dsp::cuEKF::GetQVal() {

	EKF_Update_Q<<<1,64,sizeof(double), *cuStream>>>(upQ_LPFv_vals_d, upQ_LPFv_idx_d, upQ_LPFv_avg_d, *Q.elements, Q.height, *X_k1k1.elements);
	cuBLASCheckMSt(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, F.height, Q.width, F.width,
			&constPosOne, *F.elements, F.height, *Q.elements, Q.height, &constZero, *Q_temp.elements, Q_temp.height));
	cuBLASCheckMSt(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_T, Q_temp.height, F.width, Q_temp.width,
			&constPosOne, *Q_temp.elements, Q_temp.height, *F.elements, F.height, &constZero, *Q.elements, Q.height));

	return 0;
}



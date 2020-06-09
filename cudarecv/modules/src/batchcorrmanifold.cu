

// If the grid is not correctly configured compared to the settings in main.cu,
// or if you want different spacings between points, go to:
//
//*** HACK TO LOCALLY CONFIGURE THE GRID HERE ***
//
// to strongarm the grid to the values you want.
//



#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <cuComplex.h>
#include <cufft.h>
#include <thrust/device_ptr.h>
#include <thrust/extrema.h>
#include <iostream>
#include <fstream> // Just for debug to print the grid
#include <cmath>
#include <cstring>
#include <cstdio>
#include <stdint.h>
#include "batchcorrmanifold.h"
#include "errorhandler.h"
#include <string.h>




/*
 *
 * DEVICE CONSTANTS
 *
 */

// Number of grid points is not allowed to change dynamically during the flow, so store its size as device constants
// (dynamic allocation of device memory for variable grid size could be implemented,
//  but was deemed not worth the time for initial release)
// However, support is included for dynamic changing of grid spacing (though functionality is still primitive)
__device__ __constant__ int BCM_GRID_DIMS_d[2*BCM_NUM_GRID_STATES];
__device__ __constant__ int BCM_HALF_IDX_d[2*BCM_NUM_GRID_STATES];
__device__ __constant__ int BCM_GRID_SIZE_d[2];
__device__ __constant__ int BCM_SATPOS_ARR_SIZE_d;
__device__ __constant__ int BCM_SATPOS_NUM_ARRS_d;
__device__ __constant__ double BCM_SATPOS_BATCH_TIME_d;
__device__ __constant__ int BCM_LPOWER_d;




/*
 *
 * DEVICE FUNCTIONS
 *
 */

static __device__ __host__ inline cufftDoubleComplex
BCM_ComplexScale(cufftDoubleComplex a, double s) {
	cufftDoubleComplex c;
	c.x = s * a.x;
	c.y = s * a.y;
	return c;
}




// Parallelized lat-lon calculator that returns values in radians
// Expected to be used as inline within BCM to not clutter BCM kernels
__device__ double
BCM_Dev_ECEF2LL_Rad(int idx, const double *posECEF) {
	if (idx) {
		// Compute lon
		return atan2(posECEF[1], posECEF[0]);
	}
	else {
		// Compute lat
		// NOTE: This is the "closed-form" calculation for consistency with PyGNSS (https://microem.ru/files/2012/08/GPS.G1-X-00006.pdf)
		double p = norm(2, posECEF);
		double theta = atan2(posECEF[2]*CONST_WGS84_A, p*CONST_WGS84_B);
		return atan2((posECEF[2] + pow(CONST_WGS84_EP,2) * CONST_WGS84_B * pow(sin(theta),3)),
					  (p - pow(CONST_WGS84_E,2) * CONST_WGS84_A * pow(cos(theta), 3)));
	}
}


// Given latitude and longitude, return the elements of the ENU->ECEF rotation matrix
__device__ double
BCM_Dev_Pop_R_ENU2ECEF(int idx, double *posLL) {
	switch(idx) {
		case 0:
			return -sin(posLL[1]);
		case 1:
			return -sin(posLL[0])*cos(posLL[1]);
		case 2:
			return cos(posLL[0])*cos(posLL[1]);
		case 3:
			return cos(posLL[1]);
		case 4:
			return -sin(posLL[0])*sin(posLL[1]);
		case 5:
			return	cos(posLL[0])*sin(posLL[1]);
		case 6:
			return 0.0;
		case 7:
			return cos(posLL[0]);
		case 8:
			return sin(posLL[0]);
		default:
			// TODO: tell the user that they shouldn't be here
			// (In the mean time) gum up the works
			return nan(0);
	}
}










/*
 *
 * DEVICE KERNELS
 *
 */

// Note: initializing template grids doesn't guarantee display-oriented grids!!!
//       this is basically just running linspace parallelized only once
//       so, the values produced here should be treated as dx meters in ENU (needing converted to ECEF)

/** \brief Generate the dx values of every state on the position grid
 *
 * Note: velocity grids assume that the time-domain states are symmetric about a center state!
 * Refer to VelMeas functions if you want to change this.
 *
 * \param grid				Output: The dx value of every state on the grid from the center state
 * \param gridType			Select how the grid will be constructed
 * \param gridUnitSpacing	Used in grid construction to determine spacing between states
 * \param timeGrid			Output: The unique time states on the grid -- used later when computing satellite states
 *
 */
__global__ void
BCM_InitPosGrid(dsp::utils::statePosManifold_t<double> *grid, dsp::utils::ManifoldGridTypes gridType, double *gridUnitSpacing, double *timeGrid)
{

	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int stride = blockDim.x * gridDim.x;

    int currXIdx;
    int currYIdx;
    int currZIdx;
    int currTIdx;
    int currTemp;

    bool temp;

    while (i < BCM_GRID_SIZE_d[0]) {
    	// Find the state for this i-value (assuming indexing is grid[xidx][yidx][zidx][tidx])
		currTIdx = POSMOD(i, BCM_GRID_DIMS_d[3]);
		currTemp = i / BCM_GRID_DIMS_d[3];
		currZIdx = POSMOD(currTemp, BCM_GRID_DIMS_d[2]);
		currTemp = currTemp / BCM_GRID_DIMS_d[2];
		currYIdx = POSMOD(currTemp, BCM_GRID_DIMS_d[1]);
		currXIdx = currTemp / BCM_GRID_DIMS_d[1];
		if (currXIdx > BCM_GRID_DIMS_d[0]) {
			temp = 1; // Debug
			return;
		}

    	switch(gridType) {
    		case dsp::utils::ManifoldGridTypes::Uniform:
    			// Find how far away this index is from the center point and multiply uniformly by the scale factor
    			grid[i].x = 		gridUnitSpacing[0] * (currXIdx - BCM_HALF_IDX_d[0]);
    			grid[i].y = 		gridUnitSpacing[1] * (currYIdx - BCM_HALF_IDX_d[1]);
    			grid[i].z = 		gridUnitSpacing[2] * (currZIdx - BCM_HALF_IDX_d[2]);
    			grid[i].delta_t = 	gridUnitSpacing[3] * (currTIdx - BCM_HALF_IDX_d[3]); // Note this is distance (delta_t times C)

    			// Save the time grid values during the first pass through the time states
    			if (i == currTIdx) {
    				timeGrid[i] = grid[i].delta_t;
    			}
    			break;

    		case dsp::utils::ManifoldGridTypes::ArthurBasis:
    			// Find how far away this index is from the center point and multiply uniformly by the scale factor
    			if (currXIdx < BCM_HALF_IDX_d[0]/2 || (BCM_GRID_DIMS_d[0]-currXIdx) < BCM_HALF_IDX_d[0]/2) {
    				if (currXIdx < BCM_HALF_IDX_d[0]) {
    					grid[i].x = 3 * gridUnitSpacing[0] * (currXIdx - BCM_HALF_IDX_d[0]) + gridUnitSpacing[0] * ((BCM_HALF_IDX_d[0]/2)+1) * 2;
    				}
    				else {
    					grid[i].x = 3 * gridUnitSpacing[0] * (currXIdx - BCM_HALF_IDX_d[0]) - gridUnitSpacing[0] * ((BCM_HALF_IDX_d[0]/2)+1) * 2;
    				}

    			}
    			else {
    				grid[i].x = 		gridUnitSpacing[0] * (currXIdx - BCM_HALF_IDX_d[0]);
    			}

    			if (currYIdx < BCM_HALF_IDX_d[1]/2 || (BCM_GRID_DIMS_d[1]-currYIdx) < BCM_HALF_IDX_d[1]/2) {
    				if (currYIdx < BCM_HALF_IDX_d[1]) {
    					grid[i].y = 3 * gridUnitSpacing[1] * (currYIdx - BCM_HALF_IDX_d[1]) + gridUnitSpacing[1] * ((BCM_HALF_IDX_d[1]/2)+1) * 2;
    				}
    				else {
    					grid[i].y = 3 * gridUnitSpacing[1] * (currYIdx - BCM_HALF_IDX_d[1]) - gridUnitSpacing[1] * ((BCM_HALF_IDX_d[1]/2)+1) * 2;
    				}
    			}
    			else {
    				grid[i].y = 		gridUnitSpacing[1] * (currYIdx - BCM_HALF_IDX_d[1]);
    			}

    			if (currZIdx < BCM_HALF_IDX_d[2]/2 || (BCM_GRID_DIMS_d[2]-currZIdx) < BCM_HALF_IDX_d[2]/2) {
        			if (currZIdx < BCM_HALF_IDX_d[2]) {
        				grid[i].z = 3 * gridUnitSpacing[2] * (currZIdx - BCM_HALF_IDX_d[2]) + gridUnitSpacing[2] * ((BCM_HALF_IDX_d[2]/2)+1) * 2;
        			}
        			else {
        				grid[i].z = 3 * gridUnitSpacing[2] * (currZIdx - BCM_HALF_IDX_d[2]) - gridUnitSpacing[2] * ((BCM_HALF_IDX_d[2]/2)+1) * 2;
        			}
    			}
    			else {
    				grid[i].z = 		gridUnitSpacing[2] * (currZIdx - BCM_HALF_IDX_d[2]);

    			}

    			if (currTIdx < BCM_HALF_IDX_d[3]/2 || (BCM_GRID_DIMS_d[3]-currTIdx) < BCM_HALF_IDX_d[3]/2) {
        			if (currTIdx < BCM_HALF_IDX_d[3]) {
        				grid[i].delta_t = 3 * gridUnitSpacing[3] * (currTIdx - BCM_HALF_IDX_d[3]) + gridUnitSpacing[3] * ((BCM_HALF_IDX_d[3]/2)+1) * 2;
        			}
        			else {
        				grid[i].delta_t = 3 * gridUnitSpacing[3] * (currTIdx - BCM_HALF_IDX_d[3]) - gridUnitSpacing[3] * ((BCM_HALF_IDX_d[3]/2)+1) * 2;
        			}
    			}
    			else {
    				grid[i].delta_t = 	gridUnitSpacing[3] * (currTIdx - BCM_HALF_IDX_d[3]);
    			}

    			// Save the time grid values during the first pass through the time states
    			if (i == currTIdx) {
    				timeGrid[i] = grid[i].delta_t;
    			}
    			break;
    		default:
    			// TODO: report unsupported manifold type here
    			return;
    	}

    	i += stride;

    }
}


/** \brief Generate the dx values of every state on the velocity grid
 *
 * \param grid				Output: The dx value of every state on the grid from the center state
 * \param gridType			Select how the grid will be constructed
 * \param gridUnitSpacing	Used in grid construction to determine spacing between states
 *
 */
__global__ void
BCM_InitVelGrid(dsp::utils::stateVelManifold_t<double> *grid, dsp::utils::ManifoldGridTypes gridType, double *gridUnitSpacing)
{

	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int stride = blockDim.x * gridDim.x;

    int currXIdx;
    int currYIdx;
    int currZIdx;
    int currTIdx;
    int currTemp;

    bool temp;

    while (i < BCM_GRID_SIZE_d[1]) {
    	// Find the state for this i-value (assuming indexing is grid[xidx][yidx][zidx][tidx])
		currTIdx = POSMOD(i, BCM_GRID_DIMS_d[7]);
		currTemp = i / BCM_GRID_DIMS_d[7];
		currZIdx = POSMOD(currTemp, BCM_GRID_DIMS_d[6]);
		currTemp = currTemp / BCM_GRID_DIMS_d[6];
		currYIdx = POSMOD(currTemp, BCM_GRID_DIMS_d[5]);
		currXIdx = currTemp / BCM_GRID_DIMS_d[5];
		if (currXIdx > BCM_GRID_DIMS_d[4]) {
			temp = 1; // Debug
			return;
		}

    	switch(gridType) {
    		case dsp::utils::ManifoldGridTypes::Uniform:
    			// Find how far away this index is from the center point and multiply uniformly by the scale factor
    			grid[i].x_dot = 		gridUnitSpacing[4] * (currXIdx - BCM_HALF_IDX_d[4]);
    			grid[i].y_dot = 		gridUnitSpacing[5] * (currYIdx - BCM_HALF_IDX_d[5]);
    			grid[i].z_dot = 		gridUnitSpacing[6] * (currZIdx - BCM_HALF_IDX_d[6]);
    			grid[i].delta_t_dot = 	gridUnitSpacing[7] * (currTIdx - BCM_HALF_IDX_d[7]); // Note this is distance (delta_t times C)
    			break;
    		case dsp::utils::ManifoldGridTypes::ArthurBasis:
    			// Find how far away this index is from the center point and multiply uniformly by the scale factor
    			grid[i].x_dot = 		gridUnitSpacing[4] * (currXIdx - BCM_HALF_IDX_d[4]);
    			grid[i].y_dot = 		gridUnitSpacing[5] * (currYIdx - BCM_HALF_IDX_d[5]);
    			grid[i].z_dot = 		gridUnitSpacing[6] * (currZIdx - BCM_HALF_IDX_d[6]);
    			grid[i].delta_t_dot = 	gridUnitSpacing[7] * (currTIdx - BCM_HALF_IDX_d[7]); // Note this is distance (delta_t times C)
    			break;
    		default:
    			// TODO: report unsupported manifold type here
    			return;
    	}

    	i += stride;

    }
}



/** \brief Compute the position-domain measurement as the score-weighted average state from the manifold grid using atomic adds (slow)
 *
 * Note: Expecting a 1D kernel launch
 *
 * Note: SatPos computes the relativistic corrections and factors them in to the txTime provided.
 * So, when looking up a satellite position by time, it should be looked up WITHOUT application
 * of relativistic corrections (batch SatPos effectively calculates PyGNSS's
 * satellite_clock_correction immediately followed by locate_satellite).
 * However, when computing back-calculation equations, the relativistic corrections do need to be
 * applied appropriately. These can be retrieved from the delta_t and delta_t_dot elements of
 * a given satellite state.
 *
 *
 * \param satStates 	Ptr to the state of each satellite
 * \param codeScores 	Ptr to the array of code scores (1D array of scores, with tracked PRNs adjacent)
 * \param centerPt   	Ptr to the array of the ECEF location by which the template ENU grid should be
 *        		        translated and rotated aka X_t|t-1 EKF estimate
 * \param grid	   		The template ENU grid
 * \param codeFreq   	Ptr to array of the codeFreqs of the tracked channels (used to find code scores idx)
 * \param txTime		Ptr to array of the transmission times from sat as calculated by the chan params
 * \param rxTime		The current time of the receiver's internal clock (no delta-t applied)
 * \param numChan		The number of satellites being tracked
 * \param fs			The sampling frequency of the sample set (used to find code scores idx)
 * \param numSamps		The number of samples loaded this iteration (equal to the size of one PRN of code scores)
 * \param zVal			Ptr to the array holding the measurement for this iteration
 * \param RVal			Ptr to the array holding the covariance of the measurement for this iteration
 * \param gridScores  	Ptr to array of scores for every point on the grid (for display purposes only)
 */
__global__ void
BCM_PosMeas(dsp::utils::state_t<double> *satStates, cufftDoubleComplex *codeScores,
		const double *centerPt, dsp::utils::statePosManifold_t<double> *grid,
		double *codeFreq, double *txTime, double rxTime, int numChan,
		double fs, int numSamps,
		double *zVal, double *RVal, double *gridScores)
{
	// Thread indexing variables
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int stride = blockDim.x * gridDim.x;
    unsigned int blockThreadidx = threadIdx.x;
    unsigned int blockThreadStride = blockDim.x;

    // Intermediate variables
    __shared__ double centerPtLL[2];	// lat-lon form of centerPt
    __shared__ double R_ENU2ECEF[3][3];	// Rotation matrix used by template grid so points align with ENU plots

    // Shared variables to assess the manifold (for speed) --
    __shared__ double weightedPt[4];
    __shared__ double scoresSum;

    dsp::utils::state_t<double> currPt;	// The grid point to be evaluated
    double currScore;					// The score of the current point to be evaluated


    double satPosTransmitTime;
    dsp::utils::state_t<double> currPtSatState;


    // Back-calculation values for finding scores
    double bc_rangetime;
    double bc_transmitTime;
    double bc_rc0;
    double bc_rc0_idx_base;
    double bc_rc0_idx;
    int bc_rc0_cidx;
    int bc_rc0_fidx;
    cufftDoubleComplex bc_rc0_corr;
    double cos_tau_OEDot;
    double sin_tau_OEDot;

    bool temp;
    int test=1;
    double meanPos[4]; // Debug



    // Compute the lat-lon form of centerPt and inititalize shared memory
    while (blockThreadidx < 7) {
    	switch (blockThreadidx) {
    		case 0:
    		case 1:
    			centerPtLL[blockThreadidx] = BCM_Dev_ECEF2LL_Rad(blockThreadidx, centerPt);
    			break;
    		case 2:
    		case 3:
    		case 4:
    		case 5:
    			weightedPt[blockThreadidx-2] = 0.0;
    			break;
    		case 6:
    			scoresSum = 0.0;
    			break;
    		default:
    			break;
    	}
    	blockThreadidx += blockThreadStride;
    }

    blockThreadidx = threadIdx.x;

    __syncthreads();



    // Determine rotation matrix for the template grid
    while (blockThreadidx < 9) {
    	*((double*)R_ENU2ECEF + blockThreadidx) = BCM_Dev_Pop_R_ENU2ECEF(blockThreadidx, centerPtLL); // Using pointers to treat R like 1D array while keeping 2D accessibility

    	blockThreadidx += blockThreadStride;
    }

    // Reset thread index so the next segment will start clean
    blockThreadidx = threadIdx.x;

    __syncthreads();


    // Determine scores for each candidate point
    while (i < BCM_GRID_SIZE_d[0]) {

    	// Reset the score holder for this thread
    	currScore = 0;

    	// Find the ECEF position of this candidate point
    	currPt.x = R_ENU2ECEF[0][0]*grid[i].x + R_ENU2ECEF[0][1]*grid[i].y + R_ENU2ECEF[0][2]*grid[i].z + centerPt[0];
    	currPt.y = R_ENU2ECEF[1][0]*grid[i].x + R_ENU2ECEF[1][1]*grid[i].y + R_ENU2ECEF[1][2]*grid[i].z + centerPt[1];
    	currPt.z = R_ENU2ECEF[2][0]*grid[i].x + R_ENU2ECEF[2][1]*grid[i].y + R_ENU2ECEF[2][2]*grid[i].z + centerPt[2];
    	currPt.delta_t = grid[i].delta_t + centerPt[3]; // Note that this is distance (delta_t times C)


    	// Iterate over every satellite with valid ephems
    	currScore = 0;
    	for (int currChan = 0; currChan < numChan; currChan++) {

    		// For this satellite position, compute the transmit time (TOF) to the candidate point
    		satPosTransmitTime = rxTime - (txTime[currChan] + (currPt.delta_t/(double)CONST_C)) + satStates[currChan].delta_t;

    		// Convert the satellite and coordinate positions to ECI to add
    		cos_tau_OEDot = cos(-CONST_OEDot*satPosTransmitTime);
    		sin_tau_OEDot = sin(-CONST_OEDot*satPosTransmitTime);


    		// Rotate the satellite position over the earth for this point's pvt state
    		currPtSatState.x 			= cos_tau_OEDot * satStates[currChan].x
    									- sin_tau_OEDot * satStates[currChan].y;
    		currPtSatState.y 			= sin_tau_OEDot * satStates[currChan].x
    									+ cos_tau_OEDot * satStates[currChan].y;
    		currPtSatState.z 			= satStates[currChan].z;
    		currPtSatState.delta_t 		= satStates[currChan].delta_t;


			// Find the code score index
			bc_rangetime = norm3d(currPtSatState.x - currPt.x, currPtSatState.y - currPt.y, currPtSatState.z - currPt.z) / CONST_C;
			bc_rc0 = CONST_F_CA * (satPosTransmitTime - bc_rangetime);
			bc_rc0_idx_base = (fs/codeFreq[currChan])*(-bc_rc0) + numSamps/2.0;

			// Make sure the code score index is within the range of this PRN
			// (Note: strictly less/greater so cel/flr stays within range)
			if (bc_rc0_idx_base < numSamps && bc_rc0_idx_base > 0) {
				// Also, add in the offset corresponding to this prn
				bc_rc0_idx = bc_rc0_idx_base + (numSamps*currChan);
				bc_rc0_fidx = floor(bc_rc0_idx);
				bc_rc0_cidx = floor(bc_rc0_idx+1);
			}
			else {
				//TODO: return error
				test = -4;
			}

			// Interpolate the score and add to the running total for this point
			// TODO: Do these tests need to be here? Had warp_illegal_address without them...
			cufftDoubleComplex test1 = codeScores[bc_rc0_cidx];
			cufftDoubleComplex test2 = codeScores[bc_rc0_fidx];
			bc_rc0_corr = cuCadd(BCM_ComplexScale(test1, bc_rc0_idx-bc_rc0_fidx),
								 BCM_ComplexScale(test2, bc_rc0_cidx-bc_rc0_idx));
			//bc_rc0_corr = cuCadd(BCM_ComplexScale(codeScores[bc_rc0_cidx], bc_rc0_idx-bc_rc0_fidx),
			//					 BCM_ComplexScale(codeScores[bc_rc0_fidx], bc_rc0_cidx-bc_rc0_idx));

			currScore += pow(cuCabs(bc_rc0_corr), BCM_LPOWER_d);

    	} // Done iterating over all PRNs


		// Take the score for this point and weighted-add it to the running total
    	// Note: if atomics are slow, consider a scan implementation
		atomicAdd(&weightedPt[0], currScore * currPt.x);
		atomicAdd(&weightedPt[1], currScore * currPt.y);
		atomicAdd(&weightedPt[2], currScore * currPt.z);
		atomicAdd(&weightedPt[3], currScore * currPt.delta_t);
		atomicAdd(&scoresSum, currScore);

    	// Copy score to be accessible for display (note: this is only for display, unnecessary for functionality)
    	gridScores[i] = currScore;

    	i += stride;
    }


    // Reset thread index so the next segment will start clean
    i = blockIdx.x * blockDim.x + threadIdx.x;

    __syncthreads();


    // Compute the outputs
    while (i < (32+4)) {
    	// Generate the covariance matrix as identity (for now)
    	if (i == 0 || i == 9 || i == 18 || i == 27) { RVal[i] = 1; }
    	else if (i < 32) 							{ RVal[i] = 0; }

    	// Compute the z-vals
    	else
    	{ zVal[i-32] = weightedPt[i-32] / scoresSum; }

    	i += stride;
    }
}






/** \brief Compute the velocity-domain measurement as the score-weighted average state from the manifold grid using atomic adds (slow)
 *
 * Note: Expecting a 1D kernel launch
 *
 * Note: SatPos computes the relativistic corrections and factors them in to the txTime provided.
 * So, when looking up a satellite position by time, it should be looked up WITHOUT application
 * of relativistic corrections (batch SatPos effectively calculates PyGNSS's
 * satellite_clock_correction immediately followed by locate_satellite).
 * However, when computing back-calculation equations, the relativistic corrections do need to be
 * applied appropriately. These can be retrieved from the delta_t and delta_t_dot elements of
 * a given satellite state.
 *
 *
 * \param satStates 	Ptr to the state of each satellite
 * \param carrScores 	Ptr to the array of carrier scores (1D array of scores, with tracked PRNs adjacent)
 * \param centerPt   	Ptr to the array of the ECEF location by which the template ENU grid should be
 *        		        translated and rotated aka X_t|t-1 EKF estimate
 * \param grid	   		The template ENU grid
 * \param carrFreq   	Ptr to array of the carrFreqs of the tracked channels (used to find code scores idx)
 * \param txTime		Ptr to array of the transmission times from sat as calculated by the chan params
 * \param rxTime		The current time of the receiver's internal clock (no delta-t applied)
 * \param numChan		The number of satellites being tracked
 * \param fs			The sampling frequency of the sample set (used to find code scores idx)
 * \param numfftPts		The size of one channel's worth of scores
 * \param dopplerSign	The sign convention for Doppler frequency
 * \param zVal			Ptr to the array holding the measurement for this iteration
 * \param RVal			Ptr to the array holding the covariance of the measurement for this iteration
 * \param gridScores  	Ptr to array of scores for every point on the grid (for display purposes only)
 */
__global__ void
BCM_VelMeas(dsp::utils::state_t<double> *satStates, cufftDoubleComplex *carrScores,
		const double *centerPt, dsp::utils::stateVelManifold_t<double> *grid,
		double *carrFreq, double *txTime, double rxTime, int numChan,
		double fs, int numfftPts, int *dopplerSign,
		double *zVal, double *RVal, double *gridScores)
{
	// Thread indexing variables
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int stride = blockDim.x * gridDim.x;

    // Intermediate variables
    __shared__ double rxPtLL[2];	// lat-lon form of rxPt
    __shared__ double R_ENU2ECEF[3][3];	// Rotation matrix used by template grid so points align with ENU plots

    // Shared variables to assess the manifold (for speed) --
    __shared__ double weightedPt[4];
    __shared__ double scoresSum;

    dsp::utils::stateVelManifold_t<double> currPt;	// The grid point to be evaluated
    double currScore;								// The score of the current point to be evaluated

    double satPosTransmitTime;
    dsp::utils::state_t<double> currPtSatState;

    // Back-calculation values for finding scores
    double currPtVelECI[4];
    double bc_los[3];
    double bc_range;
    double bc_losrangerate;
    double bc_fi;
    double bc_fi0;
    double bc_fi0_idx_base;
    double bc_fi0_idx;
    int bc_fi0_fidx;
    int bc_fi0_cidx;
    cufftDoubleComplex bc_fi0_fft;
    double cos_tau_OEDot;
    double sin_tau_OEDot;


    bool temp;
    double meanVel[4];
    int test = 1;


    // Note that, for this function, centerPt is a velocity point while rxPt is a position point.
    // This means that the template grid should be rotated according to rxPt so the estimates
    // still align with the ENU coordinates corresponding to the receiver position

    // Compute the lat-lon form of rxPt and initialize shared memory
    while (i < 7) {
    	switch (i) {
    		case 0:
    		case 1:
    			rxPtLL[i] = BCM_Dev_ECEF2LL_Rad(i, centerPt);
    			break;
    		case 2:
    		case 3:
    		case 4:
    		case 5:
    			weightedPt[i-2] = 0.0;
    			break;
    		case 6:
    			scoresSum = 0.0;
    			break;
    		default:
    			break;
    	}
		i += stride;
    }

    i = blockIdx.x * blockDim.x + threadIdx.x;

    __syncthreads();



    // Determine rotation matrix for the template grid
    while (i < 9) {
    	*((double*)R_ENU2ECEF + i) = BCM_Dev_Pop_R_ENU2ECEF(i, rxPtLL); // Using pointers to treat R like 1D array while keeping 2D accessibility
    	i += stride;
    }

    // Reset thread index so the next segment will start clean
    i = blockIdx.x * blockDim.x + threadIdx.x;

    __syncthreads();



    // Determine scores for each candidate point
    while (i < BCM_GRID_SIZE_d[1]) {

    	// Reset the score holder for this thread
    	currScore = 0;

    	// Find the ECEF position of this candidate point (note: this is actually a velocity state)
    	currPt.x_dot = R_ENU2ECEF[0][0]*grid[i].x_dot + R_ENU2ECEF[0][1]*grid[i].y_dot + R_ENU2ECEF[0][2]*grid[i].z_dot + centerPt[4];
    	currPt.y_dot = R_ENU2ECEF[1][0]*grid[i].x_dot + R_ENU2ECEF[1][1]*grid[i].y_dot + R_ENU2ECEF[1][2]*grid[i].z_dot + centerPt[5];
    	currPt.z_dot = R_ENU2ECEF[2][0]*grid[i].x_dot + R_ENU2ECEF[2][1]*grid[i].y_dot + R_ENU2ECEF[2][2]*grid[i].z_dot + centerPt[6];
    	currPt.delta_t_dot = grid[i].delta_t_dot + centerPt[7];


    	// Iterate over every satellite with valid ephems
    	currScore = 0;
    	for (int currChan = 0; currChan < numChan; currChan++) {

    		// For this satellite position, compute the transmit time (TOF) to the candidate point
    		satPosTransmitTime = rxTime - (txTime[currChan] + (centerPt[3]/(double)CONST_C)) + satStates[currChan].delta_t;

    		// Convert the satellite and coordinate positions to ECI to add
    		cos_tau_OEDot = cos(-CONST_OEDot*satPosTransmitTime);
    		sin_tau_OEDot = sin(-CONST_OEDot*satPosTransmitTime);


    		// Rotate the satellite position over the earth for this point's pvt state
    		currPtSatState.x 			= cos_tau_OEDot * satStates[currChan].x
    									- sin_tau_OEDot * satStates[currChan].y;
    		currPtSatState.y 			= sin_tau_OEDot * satStates[currChan].x
    									+ cos_tau_OEDot * satStates[currChan].y;
    		currPtSatState.z 			= satStates[currChan].z;
    		currPtSatState.delta_t 		= satStates[currChan].delta_t;

    		currPtSatState.x_dot		= cos_tau_OEDot * satStates[currChan].x_dot
    									- sin_tau_OEDot * satStates[currChan].y_dot
    									- CONST_OEDot * sin_tau_OEDot * satStates[currChan].x
    									- CONST_OEDot * cos_tau_OEDot * satStates[currChan].y;
    		currPtSatState.y_dot		= sin_tau_OEDot * satStates[currChan].x_dot
    									+ cos_tau_OEDot * satStates[currChan].y_dot
    									+ CONST_OEDot * cos_tau_OEDot * satStates[currChan].x
    									- CONST_OEDot * sin_tau_OEDot * satStates[currChan].y;
    		currPtSatState.z_dot		= satStates[currChan].z_dot;
    		currPtSatState.delta_t_dot	= satStates[currChan].delta_t_dot;


    		// Also need to rotate the point's velocity into the inertial frame
    		currPtVelECI[0] = currPt.x_dot - CONST_OEDot*centerPt[1];
    		currPtVelECI[1] = currPt.y_dot + CONST_OEDot*centerPt[0];
    		currPtVelECI[2] = currPt.z_dot;
    		currPtVelECI[3] = currPt.delta_t_dot;



			// Find carrier frequency (using centerPt since position states are fixed in vel manifold)
			bc_los[0] = currPtSatState.x - centerPt[0];
			bc_los[1] = currPtSatState.y - centerPt[1];
			bc_los[2] = currPtSatState.z - centerPt[2];
			bc_range = norm(3, bc_los);
			bc_losrangerate = 	((bc_los[0]/bc_range)*(currPtVelECI[0]-currPtSatState.x_dot)) +
								((bc_los[1]/bc_range)*(currPtVelECI[1]-currPtSatState.y_dot)) +
								((bc_los[2]/bc_range)*(currPtVelECI[2]-currPtSatState.z_dot));
			bc_fi = CONST_F_L1 * ((bc_losrangerate - currPtVelECI[3])/CONST_C + currPtSatState.delta_t_dot) / (*dopplerSign);


			bc_fi0 = bc_fi - carrFreq[currChan];
			// As opposed to bc_rc0, which would be negative here, based on the way fi is defined,
			// bc_fi0 should be positive here
			bc_fi0_idx_base = (numfftPts/fs)*(bc_fi0) + numfftPts/2.0;

			// Make sure the code score index is within the range of this PRN
			// (Note: strictly less/greater so cel/flr stays within range)
			if (bc_fi0_idx_base < numfftPts && bc_fi0_idx_base > 0) {
				// Also, add in the offset corresponding to this prn
				bc_fi0_idx = bc_fi0_idx_base + (numfftPts*currChan);
				bc_fi0_fidx = floor(bc_fi0_idx);
				bc_fi0_cidx = floor(bc_fi0_idx+1);
			}
			else {
				//TODO: return error
				test = -4;
			}

			// Interpolate the score and add to the running total for this point
			bc_fi0_fft = cuCadd(BCM_ComplexScale(carrScores[bc_fi0_cidx], bc_fi0_idx-bc_fi0_fidx),
								BCM_ComplexScale(carrScores[bc_fi0_fidx], bc_fi0_cidx-bc_fi0_idx));
			currScore += pow(cuCabs(bc_fi0_fft), BCM_LPOWER_d);

    	} // Done iterating over all channels

		// Take the score for this point and weighted-add it to the running total
		atomicAdd(&weightedPt[0], currScore * currPt.x_dot);
		atomicAdd(&weightedPt[1], currScore * currPt.y_dot);
		atomicAdd(&weightedPt[2], currScore * currPt.z_dot);
		atomicAdd(&weightedPt[3], currScore * currPt.delta_t_dot);
		atomicAdd(&scoresSum, currScore);

    	// Copy score to be accessible for display (note: this is only for display, unnecessary for functionality)
    	gridScores[i] = currScore;

    	i += stride;
    }


    // Reset thread index so the next segment will start clean
    i = blockIdx.x * blockDim.x + threadIdx.x;

    __syncthreads();


    // Compute the outputs
    while (i < (32+4)) {
    	// Generate the covariance matrix as identity (for now) with +4 offset so the identity continues on row 5
    	if (i == (0+4) || i == (9+4) || i == (18+4) || i == (27+4)) { RVal[i+32] = 1; }
    	else if (i < 32) 											{ RVal[i+32] = 0; }

    	// Compute the z-vals
    	else
    	{ zVal[i-32+4] = weightedPt[i-32] / scoresSum; }

    	i += stride;
    }

}




/** \brief Compute the position-domain measurement as the score-weighted average state from the manifold grid using parallel reductions
 *
 * Note: Expecting a 1D kernel launch
 *
 * Note: SatPos computes the relativistic corrections and factors them in to the txTime provided.
 * So, when looking up a satellite position by time, it should be looked up WITHOUT application
 * of relativistic corrections (batch SatPos effectively calculates PyGNSS's
 * satellite_clock_correction immediately followed by locate_satellite).
 * However, when computing back-calculation equations, the relativistic corrections do need to be
 * applied appropriately. These can be retrieved from the delta_t and delta_t_dot elements of
 * a given satellite state.
 *
 *
 * \param satStates 		Ptr to the state of each satellite
 * \param codeScores 		Ptr to the array of code scores (1D array of scores, with tracked PRNs adjacent)
 * \param centerPt   		Ptr to the array of the ECEF location by which the template ENU grid should be
 *        		        	translated and rotated aka X_t|t-1 EKF estimate
 * \param grid	   			The template ENU grid
 * \param enu2ecefMat		Ptr to a rotation matrix from ENU (grid domain) to ECEF (satellite domain)
 * \param codeFreq   		Ptr to array of the codeFreqs of the tracked channels (used to find code scores idx)
 * \param txTime			Ptr to array of the transmission times from sat as calculated by the chan params
 * \param rxTime			The current time of the receiver's internal clock (no delta-t applied)
 * \param numChan			The number of satellites being tracked
 * \param fs				The sampling frequency of the sample set (used to find code scores idx)
 * \param numSamps			The number of samples loaded this iteration (equal to the size of one PRN of code scores)
 * \param weightStateArr	Output: The measurement from this sampleset
 */
__global__ void
BCM_PosMeasReduction(const dsp::utils::state_t<double> *satStates, const cufftDoubleComplex *codeScores,
		const double *centerPt, const dsp::utils::statePosManifold_t<double> *grid, const double *enu2ecefMat,
		const double *codeFreq, const double *txTime, const double rxTime, const int numChan,
		const double fs, const int numSamps,
		dsp::utils::weightState_t<double> *weightStateArr)
{
	// Thread indexing variables
	unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int stride = blockDim.x * gridDim.x;

    // Shared memory for the reduction sum
    __shared__ dsp::utils::weightState_t<double> sdataPosMani[64];
    dsp::utils::weightState_t<double> mySum = {};

    dsp::utils::state_t<double> currPt;	// The grid point to be evaluated
    double currScore;					// The score of the current point to be evaluated

    double satPosTransmitTime;
    dsp::utils::state_t<double> currPtSatState;

    // Back-calculation values for finding scores
    double bc_rangetime;
    double bc_rc0;
    double bc_rc0_idx_base;
    double bc_rc0_idx;
    int bc_rc0_cidx;
    int bc_rc0_fidx;
    cufftDoubleComplex bc_rc0_corr;

    int currTIdx;

    int test = 1;


    // Determine scores for each candidate point
    while (i < BCM_GRID_SIZE_d[0]) {

    	// Reset the score holder for this thread
    	currScore = 0;

    	// Find the ECEF position of this candidate point
    	currPt.x = enu2ecefMat[0]*grid[i].x + enu2ecefMat[1]*grid[i].y + enu2ecefMat[2]*grid[i].z + centerPt[0];
    	currPt.y = enu2ecefMat[3]*grid[i].x + enu2ecefMat[4]*grid[i].y + enu2ecefMat[5]*grid[i].z + centerPt[1];
    	currPt.z = enu2ecefMat[6]*grid[i].x + enu2ecefMat[7]*grid[i].y + enu2ecefMat[8]*grid[i].z + centerPt[2];
    	currPt.delta_t = grid[i].delta_t + centerPt[3]; // Note that this is distance (delta_t times C)

    	// Determine what time state we're looking at
		currTIdx = POSMOD(i, BCM_GRID_DIMS_d[3]);


    	// Iterate over every satellite with valid ephems
    	currScore = 0;
    	for (int currChan = 0; currChan < numChan; currChan++) {

    		// Get the pre-calculated satellite state
    		currPtSatState = satStates[currChan*BCM_GRID_DIMS_d[3] + currTIdx];

    		// For this satellite position, compute the transmit time (TOF) to the candidate point
    		satPosTransmitTime = rxTime - (txTime[currChan] + (currPt.delta_t/(double)CONST_C)) + currPtSatState.delta_t;


			// Find the code score index
			bc_rangetime = norm3d(currPtSatState.x - currPt.x, currPtSatState.y - currPt.y, currPtSatState.z - currPt.z) / CONST_C;
			bc_rc0 = CONST_F_CA * (satPosTransmitTime - bc_rangetime);
			bc_rc0_idx_base = (fs/codeFreq[currChan])*(-bc_rc0) + numSamps/2.0;

			// Make sure the code score index is within the range of this PRN
			// (Note: strictly less/greater so cel/flr stays within range)
			if (bc_rc0_idx_base < numSamps && bc_rc0_idx_base > 0) {
				// Also, add in the offset corresponding to this prn
				bc_rc0_idx = bc_rc0_idx_base + (numSamps*currChan);
				bc_rc0_fidx = floor(bc_rc0_idx);
				bc_rc0_cidx = floor(bc_rc0_idx+1);
			}
			else {
				//TODO: return error
				test = -4;
			}

			// Interpolate the score and add to the running total for this point
			// TODO: Do these tests need to be here? Had warp_illegal_address without them...
			cufftDoubleComplex test1 = codeScores[bc_rc0_cidx];
			cufftDoubleComplex test2 = codeScores[bc_rc0_fidx];
			bc_rc0_corr = cuCadd(BCM_ComplexScale(test1, bc_rc0_idx-bc_rc0_fidx),
								 BCM_ComplexScale(test2, bc_rc0_cidx-bc_rc0_idx));
			//bc_rc0_corr = cuCadd(BCM_ComplexScale(codeScores[bc_rc0_cidx], bc_rc0_idx-bc_rc0_fidx),
			//					 BCM_ComplexScale(codeScores[bc_rc0_fidx], bc_rc0_cidx-bc_rc0_idx));

			currScore += pow(cuCabs(bc_rc0_corr), BCM_LPOWER_d);

    	} // Done iterating over all PRNs

    	// Add this score to the weighted point accumulation for this thread
    	mySum.a += currScore*currPt.x;
    	mySum.b += currScore*currPt.y;
    	mySum.c += currScore*currPt.z;
    	mySum.d += currScore*currPt.delta_t;
    	mySum.score += currScore;

    	i += stride;
    }

    // Copy the thread's local sum into shared memory for the reduction sum
    sdataPosMani[tid] = mySum;


    // Reset thread index so the next segment will start clean
    i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i == 0) {
    	test = 123;
    }


    __syncthreads();



    // Run a reduction sum for all threads in this block
    // Note: doing the full unroll reduction sum not __shfl_down_sync
    // because __shfl_down_sync only works on basic datatypes and not structs.

	// do reduction in shared mem
	if ((blockDim.x >= 512) && (tid < 256))
	{
		mySum.a += sdataPosMani[tid + 256].a;
		mySum.b += sdataPosMani[tid + 256].b;
		mySum.c += sdataPosMani[tid + 256].c;
		mySum.d += sdataPosMani[tid + 256].d;
		mySum.score += sdataPosMani[tid + 256].score;

		sdataPosMani[tid] = mySum;
	}

	__syncthreads();

	if ((blockDim.x >= 256) &&(tid < 128))
	{
		mySum.a += sdataPosMani[tid + 128].a;
		mySum.b += sdataPosMani[tid + 128].b;
		mySum.c += sdataPosMani[tid + 128].c;
		mySum.d += sdataPosMani[tid + 128].d;
		mySum.score += sdataPosMani[tid + 128].score;

		sdataPosMani[tid] = mySum;
	}

	 __syncthreads();

	if ((blockDim.x >= 128) && (tid <  64))
	{
		mySum.a += sdataPosMani[tid + 64].a;
		mySum.b += sdataPosMani[tid + 64].b;
		mySum.c += sdataPosMani[tid + 64].c;
		mySum.d += sdataPosMani[tid + 64].d;
		mySum.score += sdataPosMani[tid + 64].score;

		sdataPosMani[tid] = mySum;
	}

	__syncthreads();

	if ((blockDim.x >= 64) && (tid <  32))
	{
		mySum.a += sdataPosMani[tid + 32].a;
		mySum.b += sdataPosMani[tid + 32].b;
		mySum.c += sdataPosMani[tid + 32].c;
		mySum.d += sdataPosMani[tid + 32].d;
		mySum.score += sdataPosMani[tid + 32].score;

		sdataPosMani[tid] = mySum;
	}

	__syncthreads();

	if ((blockDim.x >= 32) && (tid <  16))
	{
		mySum.a += sdataPosMani[tid + 16].a;
		mySum.b += sdataPosMani[tid + 16].b;
		mySum.c += sdataPosMani[tid + 16].c;
		mySum.d += sdataPosMani[tid + 16].d;
		mySum.score += sdataPosMani[tid + 16].score;

		sdataPosMani[tid] = mySum;
	}

	__syncthreads();

	if ((blockDim.x >= 16) && (tid <  8))
	{
		mySum.a += sdataPosMani[tid + 8].a;
		mySum.b += sdataPosMani[tid + 8].b;
		mySum.c += sdataPosMani[tid + 8].c;
		mySum.d += sdataPosMani[tid + 8].d;
		mySum.score += sdataPosMani[tid + 8].score;

		sdataPosMani[tid] = mySum;
	}

	__syncthreads();

	if ((blockDim.x >= 8) && (tid <  4))
	{
		mySum.a += sdataPosMani[tid + 4].a;
		mySum.b += sdataPosMani[tid + 4].b;
		mySum.c += sdataPosMani[tid + 4].c;
		mySum.d += sdataPosMani[tid + 4].d;
		mySum.score += sdataPosMani[tid + 4].score;

		sdataPosMani[tid] = mySum;
	}

	__syncthreads();

	if ((blockDim.x >= 4) && (tid <  2))
	{
		mySum.a += sdataPosMani[tid + 2].a;
		mySum.b += sdataPosMani[tid + 2].b;
		mySum.c += sdataPosMani[tid + 2].c;
		mySum.d += sdataPosMani[tid + 2].d;
		mySum.score += sdataPosMani[tid + 2].score;

		sdataPosMani[tid] = mySum;
	}

	__syncthreads();

	if ((blockDim.x >= 2) && (tid <  1))
	{
		mySum.a += sdataPosMani[tid + 1].a;
		mySum.b += sdataPosMani[tid + 1].b;
		mySum.c += sdataPosMani[tid + 1].c;
		mySum.d += sdataPosMani[tid + 1].d;
		mySum.score += sdataPosMani[tid + 1].score;

		// Write result to output array
		weightStateArr[blockIdx.x] = mySum;
	}
}




/** \brief Compute the velocity-domain measurement as the score-weighted average state from the manifold grid using parallel reductions
 *
 * Note: Expecting a 1D kernel launch
 *
 * Note: SatPos computes the relativistic corrections and factors them in to the txTime provided.
 * So, when looking up a satellite position by time, it should be looked up WITHOUT application
 * of relativistic corrections (batch SatPos effectively calculates PyGNSS's
 * satellite_clock_correction immediately followed by locate_satellite).
 * However, when computing back-calculation equations, the relativistic corrections do need to be
 * applied appropriately. These can be retrieved from the delta_t and delta_t_dot elements of
 * a given satellite state.
 *
 *
 * \param satStates 		Ptr to the state of each satellite
 * \param carrScores 		Ptr to the array of carrier scores (1D array of scores, with tracked PRNs adjacent)
 * \param centerPt   		Ptr to the array of the ECEF location by which the template ENU grid should be
 *        		        	translated and rotated aka X_t|t-1 EKF estimate
 * \param grid	   			The template ENU grid
 * \param enu2ecefMat		Ptr to a rotation matrix from ENU (grid domain) to ECEF (satellite domain)
 * \param carrFreq   		Ptr to array of the carrFreqs of the tracked channels (used to find code scores idx)
 * \param txTime			Ptr to array of the transmission times from sat as calculated by the chan params
 * \param rxTime			The current time of the receiver's internal clock (no delta-t applied)
 * \param numChan			The number of satellites being tracked
 * \param fs				The sampling frequency of the sample set (used to find code scores idx)
 * \param numfftPts			The size of one channel's worth of scores
 * \param dopplerSign		The sign convention for Doppler frequency
 * \param weightStateArr	Output: The measurement from this sampleset
 *
 */
__global__ void
BCM_VelMeasReduction(const dsp::utils::state_t<double> *satStates, const cufftDoubleComplex *carrScores,
		const double *centerPt, const dsp::utils::stateVelManifold_t<double> *grid, const double *enu2ecefMat,
		const double *carrFreq, const double *txTime, const double rxTime, const int numChan,
		const double fs, const int numfftPts, const int *dopplerSign,
		dsp::utils::weightState_t<double> *weightStateArr)
{
	// Thread indexing variables
	unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int stride = blockDim.x * gridDim.x;

    // Shared memory for the reduction sum
    __shared__ dsp::utils::weightState_t<double> sdataVelMani[64];
    dsp::utils::weightState_t<double> mySum = {};

    dsp::utils::stateVelManifold_t<double> currPt;	// The grid point to be evaluated
    double currScore;								// The score of the current point to be evaluated

    dsp::utils::state_t<double> currPtSatState;

    // Back-calculation values for finding scores
    double currPtVelECI[4];
    double bc_los[3];
    double bc_range;
    double bc_losrangerate;
    double bc_fi;
    double bc_fi0;
    double bc_fi0_idx_base;
    double bc_fi0_idx;
    int bc_fi0_fidx;
    int bc_fi0_cidx;
    cufftDoubleComplex bc_fi0_fft;


    int test = 0;
    dsp::utils::weightState_t<double> testState = {};



    // Determine scores for each candidate point
    while (i < BCM_GRID_SIZE_d[1]) {

    	// Reset the score holder for this thread
    	currScore = 0;

    	// Find the ECEF position of this candidate point (note: this is actually a velocity state)
    	currPt.x_dot = enu2ecefMat[0]*grid[i].x_dot + enu2ecefMat[1]*grid[i].y_dot + enu2ecefMat[2]*grid[i].z_dot + centerPt[4];
    	currPt.y_dot = enu2ecefMat[3]*grid[i].x_dot + enu2ecefMat[4]*grid[i].y_dot + enu2ecefMat[5]*grid[i].z_dot + centerPt[5];
    	currPt.z_dot = enu2ecefMat[6]*grid[i].x_dot + enu2ecefMat[7]*grid[i].y_dot + enu2ecefMat[8]*grid[i].z_dot + centerPt[6];
    	currPt.delta_t_dot = grid[i].delta_t_dot + centerPt[7];


    	// Iterate over every satellite with valid ephems
    	currScore = 0;
    	for (int currChan = 0; currChan < numChan; currChan++) {

    		// ***** Assuming symmetric time grid, take the middle point in the set for this channel
    		currPtSatState = satStates[currChan*BCM_GRID_DIMS_d[3] + (BCM_GRID_DIMS_d[3]/2)];


    		// Also need to rotate the point's velocity into the inertial frame
    		currPtVelECI[0] = currPt.x_dot - CONST_OEDot*centerPt[1];
    		currPtVelECI[1] = currPt.y_dot + CONST_OEDot*centerPt[0];
    		currPtVelECI[2] = currPt.z_dot;
    		currPtVelECI[3] = currPt.delta_t_dot;



			// Find carrier frequency (using centerPt since position states are fixed in vel manifold)
			bc_los[0] = currPtSatState.x - centerPt[0];
			bc_los[1] = currPtSatState.y - centerPt[1];
			bc_los[2] = currPtSatState.z - centerPt[2];
			bc_range = norm(3, bc_los); // Might as well find range this way since the unit vector is needed.
			bc_losrangerate = 	((bc_los[0]/bc_range)*(currPtVelECI[0]-currPtSatState.x_dot)) +
								((bc_los[1]/bc_range)*(currPtVelECI[1]-currPtSatState.y_dot)) +
								((bc_los[2]/bc_range)*(currPtVelECI[2]-currPtSatState.z_dot));
			bc_fi = CONST_F_L1 * ((bc_losrangerate - currPtVelECI[3])/CONST_C + currPtSatState.delta_t_dot) / (*dopplerSign);


			bc_fi0 = bc_fi - carrFreq[currChan];
			// As opposed to bc_rc0, which would be negative here, based on the way fi is defined,
			// bc_fi0 should be positive here
			bc_fi0_idx_base = (numfftPts/fs)*(bc_fi0) + numfftPts/2.0;

			// Make sure the code score index is within the range of this PRN
			// (Note: strictly less/greater so cel/flr stays within range)
			if (bc_fi0_idx_base < numfftPts && bc_fi0_idx_base > 0) {
				// Also, add in the offset corresponding to this prn
				bc_fi0_idx = bc_fi0_idx_base + (numfftPts*currChan);
				bc_fi0_fidx = floor(bc_fi0_idx);
				bc_fi0_cidx = floor(bc_fi0_idx+1);
			}
			else {
				//TODO: return error
				test = -4;
			}

			// Interpolate the score and add to the running total for this point
			bc_fi0_fft = cuCadd(BCM_ComplexScale(carrScores[bc_fi0_cidx], bc_fi0_idx-bc_fi0_fidx),
								BCM_ComplexScale(carrScores[bc_fi0_fidx], bc_fi0_cidx-bc_fi0_idx));
			currScore += pow(cuCabs(bc_fi0_fft), BCM_LPOWER_d);


    	} // Done iterating over all channels

    	// Add this score to the weighted point accumulation for this thread
    	mySum.a += currScore*currPt.x_dot;
    	mySum.b += currScore*currPt.y_dot;
    	mySum.c += currScore*currPt.z_dot;
    	mySum.d += currScore*currPt.delta_t_dot;
    	mySum.score += currScore;

    	// Copy score to be accessible for display (note: this is only for display, unnecessary for functionality)
    	//gridScores[i] = currScore;

    	i += stride;
    }


    // Copy the thread's local sum into shared memory for the reduction sum
    sdataVelMani[tid] = mySum;


    // Reset thread index so the next segment will start clean
    i = blockIdx.x * blockDim.x + threadIdx.x;


    __syncthreads();



    // Run a reduction sum for all threads in this block

	// do reduction in shared mem
	if ((blockDim.x >= 512) && (tid < 256))
	{
		if (tid == 0) {
			testState = sdataVelMani[tid + 256];
			test = 29;
		}

		mySum.a += sdataVelMani[tid + 256].a;
		mySum.b += sdataVelMani[tid + 256].b;
		mySum.c += sdataVelMani[tid + 256].c;
		mySum.d += sdataVelMani[tid + 256].d;
		mySum.score += sdataVelMani[tid + 256].score;

		sdataVelMani[tid] = mySum;
	}

	__syncthreads();

	if ((blockDim.x >= 256) &&(tid < 128))
	{
		mySum.a += sdataVelMani[tid + 128].a;
		mySum.b += sdataVelMani[tid + 128].b;
		mySum.c += sdataVelMani[tid + 128].c;
		mySum.d += sdataVelMani[tid + 128].d;
		mySum.score += sdataVelMani[tid + 128].score;

		sdataVelMani[tid] = mySum;
	}

	 __syncthreads();

	if ((blockDim.x >= 128) && (tid <  64))
	{
		mySum.a += sdataVelMani[tid + 64].a;
		mySum.b += sdataVelMani[tid + 64].b;
		mySum.c += sdataVelMani[tid + 64].c;
		mySum.d += sdataVelMani[tid + 64].d;
		mySum.score += sdataVelMani[tid + 64].score;

		sdataVelMani[tid] = mySum;
	}

	__syncthreads();

	if ((blockDim.x >= 64) && (tid <  32))
	{
		mySum.a += sdataVelMani[tid + 32].a;
		mySum.b += sdataVelMani[tid + 32].b;
		mySum.c += sdataVelMani[tid + 32].c;
		mySum.d += sdataVelMani[tid + 32].d;
		mySum.score += sdataVelMani[tid + 32].score;

		sdataVelMani[tid] = mySum;
	}

	__syncthreads();

	if ((blockDim.x >= 32) && (tid <  16))
	{
		mySum.a += sdataVelMani[tid + 16].a;
		mySum.b += sdataVelMani[tid + 16].b;
		mySum.c += sdataVelMani[tid + 16].c;
		mySum.d += sdataVelMani[tid + 16].d;
		mySum.score += sdataVelMani[tid + 16].score;

		sdataVelMani[tid] = mySum;
	}

	__syncthreads();

	if ((blockDim.x >= 16) && (tid <  8))
	{
		mySum.a += sdataVelMani[tid + 8].a;
		mySum.b += sdataVelMani[tid + 8].b;
		mySum.c += sdataVelMani[tid + 8].c;
		mySum.d += sdataVelMani[tid + 8].d;
		mySum.score += sdataVelMani[tid + 8].score;

		sdataVelMani[tid] = mySum;
	}

	__syncthreads();

	if ((blockDim.x >= 8) && (tid <  4))
	{
		mySum.a += sdataVelMani[tid + 4].a;
		mySum.b += sdataVelMani[tid + 4].b;
		mySum.c += sdataVelMani[tid + 4].c;
		mySum.d += sdataVelMani[tid + 4].d;
		mySum.score += sdataVelMani[tid + 4].score;

		sdataVelMani[tid] = mySum;
	}

	__syncthreads();

	if ((blockDim.x >= 4) && (tid <  2))
	{
		mySum.a += sdataVelMani[tid + 2].a;
		mySum.b += sdataVelMani[tid + 2].b;
		mySum.c += sdataVelMani[tid + 2].c;
		mySum.d += sdataVelMani[tid + 2].d;
		mySum.score += sdataVelMani[tid + 2].score;

		sdataVelMani[tid] = mySum;
	}

	__syncthreads();

	if ((blockDim.x >= 2) && (tid <  1))
	{
		mySum.a += sdataVelMani[tid + 1].a;
		mySum.b += sdataVelMani[tid + 1].b;
		mySum.c += sdataVelMani[tid + 1].c;
		mySum.d += sdataVelMani[tid + 1].d;
		mySum.score += sdataVelMani[tid + 1].score;

		// Write result to global memory
		weightStateArr[blockIdx.x] = mySum;
	}


}







/** \brief Compute the position-domain measurement as the score-weighted average state from the manifold grid using parallel reductions
 *
 * Note: Expecting a 1D kernel launch
 *
 * \param posWeightStates	The position-domain states, weighted by their scores
 * \param numElemsToSum		The number of elements in posWeightStates
 * \param zVal				Ptr to the array holding the measurement for this iteration
 * \param RVal				Ptr to the array holding the covariance of the measurement for this iteration
 *
 */
__global__ void
BCM_ReduceAndPosMeas(dsp::utils::weightState_t<double> *posWeightStates, int numElemsToSum, double *zVal, double *RVal) {


	// Thread indexing variables
	unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int stride = blockDim.x * gridDim.x;

    int test = 0;


    // Shared memory for the reduction sum
    __shared__ dsp::utils::weightState_t<double> sdataPos[32];
    dsp::utils::weightState_t<double> mySumPos = posWeightStates[i];

	sdataPos[tid] = mySumPos;

    if (i == 0) {
    	test = 123;
    }

	if (i < numElemsToSum) {
		if ((blockDim.x >= 256) && (tid < 128))
		{
			mySumPos.a += sdataPos[tid + 128].a;
			mySumPos.b += sdataPos[tid + 128].b;
			mySumPos.c += sdataPos[tid + 128].c;
			mySumPos.d += sdataPos[tid + 128].d;
			mySumPos.score += sdataPos[tid + 128].score;

			sdataPos[tid] = mySumPos;
		}

		 __syncthreads();

		if ((blockDim.x >= 128) && (tid <  64))
		{
			mySumPos.a += sdataPos[tid + 64].a;
			mySumPos.b += sdataPos[tid + 64].b;
			mySumPos.c += sdataPos[tid + 64].c;
			mySumPos.d += sdataPos[tid + 64].d;
			mySumPos.score += sdataPos[tid + 64].score;

			sdataPos[tid] = mySumPos;
		}

		__syncthreads();

		if ((blockDim.x >= 64) && (tid <  32))
		{
			mySumPos.a += sdataPos[tid + 32].a;
			mySumPos.b += sdataPos[tid + 32].b;
			mySumPos.c += sdataPos[tid + 32].c;
			mySumPos.d += sdataPos[tid + 32].d;
			mySumPos.score += sdataPos[tid + 32].score;

			sdataPos[tid] = mySumPos;
		}

		__syncthreads();

		if ((blockDim.x >= 32) && (tid <  16))
		{
			mySumPos.a += sdataPos[tid + 16].a;
			mySumPos.b += sdataPos[tid + 16].b;
			mySumPos.c += sdataPos[tid + 16].c;
			mySumPos.d += sdataPos[tid + 16].d;
			mySumPos.score += sdataPos[tid + 16].score;

			sdataPos[tid] = mySumPos;
		}

		__syncthreads();

		if ((blockDim.x >= 16) && (tid <  8))
		{
			mySumPos.a += sdataPos[tid + 8].a;
			mySumPos.b += sdataPos[tid + 8].b;
			mySumPos.c += sdataPos[tid + 8].c;
			mySumPos.d += sdataPos[tid + 8].d;
			mySumPos.score += sdataPos[tid + 8].score;

			sdataPos[tid] = mySumPos;
		}

		__syncthreads();

		if ((blockDim.x >= 8) && (tid <  4))
		{
			mySumPos.a += sdataPos[tid + 4].a;
			mySumPos.b += sdataPos[tid + 4].b;
			mySumPos.c += sdataPos[tid + 4].c;
			mySumPos.d += sdataPos[tid + 4].d;
			mySumPos.score += sdataPos[tid + 4].score;

			sdataPos[tid] = mySumPos;
		}

		__syncthreads();

		if ((blockDim.x >= 4) && (tid <  2))
		{
			mySumPos.a += sdataPos[tid + 2].a;
			mySumPos.b += sdataPos[tid + 2].b;
			mySumPos.c += sdataPos[tid + 2].c;
			mySumPos.d += sdataPos[tid + 2].d;
			mySumPos.score += sdataPos[tid + 2].score;

			sdataPos[tid] = mySumPos;
		}

		__syncthreads();

		if ((blockDim.x >= 2) && (tid <  1))
		{
			mySumPos.a += sdataPos[tid + 1].a;
			mySumPos.b += sdataPos[tid + 1].b;
			mySumPos.c += sdataPos[tid + 1].c;
			mySumPos.d += sdataPos[tid + 1].d;
			mySumPos.score += sdataPos[tid + 1].score;

			sdataPos[tid] = mySumPos;
		}
	}

	__syncthreads();



	// Generate the measurement if you're the 0-indexed thread
	if (tid == 0) {
		zVal[0] = mySumPos.a / mySumPos.score;
		zVal[1] = mySumPos.b / mySumPos.score;
		zVal[2] = mySumPos.c / mySumPos.score;
		zVal[3] = mySumPos.d / mySumPos.score;
	}

	while (i < 32) {
		// Generate the covariance matrix as identity (for now)
		if (i == 0 || i == 9 || i == 18 || i == 27) { RVal[i] = 1; }
		else if (i < 32) 							{ RVal[i] = 0; }
		i += stride;
	}

}



/** \brief Compute the velocity-domain measurement as the score-weighted average state from the manifold grid using parallel reductions
 *
 * Note: Expecting a 1D kernel launch
 *
 * \param posWeightStates	The position-domain states, weighted by their scores
 * \param numElemsToSum		The number of elements in posWeightStates
 * \param zVal				Ptr to the array holding the measurement for this iteration
 * \param RVal				Ptr to the array holding the covariance of the measurement for this iteration
 *
 */
__global__ void
BCM_ReduceAndVelMeas(dsp::utils::weightState_t<double> *velWeightStates, int numElemsToSum, double *zVal, double *RVal) {


	// Thread indexing variables
	unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int stride = blockDim.x * gridDim.x;

    int test = 0;


    // Shared memory for the reduction sum
    __shared__ dsp::utils::weightState_t<double> sdataVel[32];
    dsp::utils::weightState_t<double> mySumVel = velWeightStates[i];

	sdataVel[tid] = mySumVel;


    if (i == 0) {
    	test = 123;
    }


    if (i < numElemsToSum) {
		if ((blockDim.x >= 256) && (tid < 128))
		{
			mySumVel.a += sdataVel[tid + 128].a;
			mySumVel.b += sdataVel[tid + 128].b;
			mySumVel.c += sdataVel[tid + 128].c;
			mySumVel.d += sdataVel[tid + 128].d;
			mySumVel.score += sdataVel[tid + 128].score;

			sdataVel[tid] = mySumVel;
		}

		 __syncthreads();

		if ((blockDim.x >= 128) && (tid <  64))
		{
			mySumVel.a += sdataVel[tid + 64].a;
			mySumVel.b += sdataVel[tid + 64].b;
			mySumVel.c += sdataVel[tid + 64].c;
			mySumVel.d += sdataVel[tid + 64].d;
			mySumVel.score += sdataVel[tid + 64].score;

			sdataVel[tid] = mySumVel;
		}

		__syncthreads();

		if ((blockDim.x >= 64) && (tid <  32))
		{
			mySumVel.a += sdataVel[tid + 32].a;
			mySumVel.b += sdataVel[tid + 32].b;
			mySumVel.c += sdataVel[tid + 32].c;
			mySumVel.d += sdataVel[tid + 32].d;
			mySumVel.score += sdataVel[tid + 32].score;

			sdataVel[tid] = mySumVel;
		}

		__syncthreads();

		if ((blockDim.x >= 32) && (tid <  16))
		{
			mySumVel.a += sdataVel[tid + 16].a;
			mySumVel.b += sdataVel[tid + 16].b;
			mySumVel.c += sdataVel[tid + 16].c;
			mySumVel.d += sdataVel[tid + 16].d;
			mySumVel.score += sdataVel[tid + 16].score;

			sdataVel[tid] = mySumVel;
		}

		__syncthreads();

		if ((blockDim.x >= 16) && (tid <  8))
		{
			mySumVel.a += sdataVel[tid + 8].a;
			mySumVel.b += sdataVel[tid + 8].b;
			mySumVel.c += sdataVel[tid + 8].c;
			mySumVel.d += sdataVel[tid + 8].d;
			mySumVel.score += sdataVel[tid + 8].score;

			sdataVel[tid] = mySumVel;
		}

		__syncthreads();

		if ((blockDim.x >= 8) && (tid <  4))
		{
			mySumVel.a += sdataVel[tid + 4].a;
			mySumVel.b += sdataVel[tid + 4].b;
			mySumVel.c += sdataVel[tid + 4].c;
			mySumVel.d += sdataVel[tid + 4].d;
			mySumVel.score += sdataVel[tid + 4].score;

			sdataVel[tid] = mySumVel;
		}

		__syncthreads();

		if ((blockDim.x >= 4) && (tid <  2))
		{
			mySumVel.a += sdataVel[tid + 2].a;
			mySumVel.b += sdataVel[tid + 2].b;
			mySumVel.c += sdataVel[tid + 2].c;
			mySumVel.d += sdataVel[tid + 2].d;
			mySumVel.score += sdataVel[tid + 2].score;

			sdataVel[tid] = mySumVel;
		}

		__syncthreads();

		if ((blockDim.x >= 2) && (tid <  1))
		{
			mySumVel.a += sdataVel[tid + 1].a;
			mySumVel.b += sdataVel[tid + 1].b;
			mySumVel.c += sdataVel[tid + 1].c;
			mySumVel.d += sdataVel[tid + 1].d;
			mySumVel.score += sdataVel[tid + 1].score;

			sdataVel[tid] = mySumVel;
		}
    }

	__syncthreads();



	// Generate the measurement if you're the 0-indexed thread
	if (tid == 0) {
		zVal[4] = mySumVel.a / mySumVel.score;
		zVal[5] = mySumVel.b / mySumVel.score;
		zVal[6] = mySumVel.c / mySumVel.score;
		zVal[7] = mySumVel.d / mySumVel.score;
	}

	while (i < 32) {
		// Generate the covariance matrix as identity (for now)
		if (i == (0+4) || i == (9+4) || i == (18+4) || i == (27+4)) { RVal[i+32] = 1; }
		else if (i < 32) 									  		{ RVal[i+32] = 0; }
		i += stride;
	}
}







/** \brief Compute the position-domain measurement as the highest score state from the manifold grid
 *
 * Note: Expecting a 1D kernel launch
 *
 * Note: SatPos computes the relativistic corrections and factors them in to the txTime provided.
 * So, when looking up a satellite position by time, it should be looked up WITHOUT application
 * of relativistic corrections (batch SatPos effectively calculates PyGNSS's
 * satellite_clock_correction immediately followed by locate_satellite).
 * However, when computing back-calculation equations, the relativistic corrections do need to be
 * applied appropriately. These can be retrieved from the delta_t and delta_t_dot elements of
 * a given satellite state.
 *
 *
 * \param satStates 		Ptr to the state of each satellite
 * \param codeScores 		Ptr to the array of code scores (1D array of scores, with tracked PRNs adjacent)
 * \param centerPt   		Ptr to the array of the ECEF location by which the template ENU grid should be
 *        		        	translated and rotated aka X_t|t-1 EKF estimate
 * \param grid	   			The template ENU grid
 * \param enu2ecefMat		Ptr to a rotation matrix from ENU (grid domain) to ECEF (satellite domain)
 * \param codeFreq   		Ptr to array of the codeFreqs of the tracked channels (used to find code scores idx)
 * \param cpRefTOW			The TOW at the reference code period for each channel
 * \param cpElapsedEnd		The number of completed code periods at the end of this sampleset for each channel
 * \param cpRef				The reference code period for each channel
 * \param codePhase			The current code phase estimate
 * \param txTime			Ptr to array of the transmission times from sat as calculated by the chan params
 * \param rxTime			The current time of the receiver's internal clock (no delta-t applied)
 * \param numChan			The number of satellites being tracked
 * \param fs				The sampling frequency of the sample set (used to find code scores idx)
 * \param numSamps			The size of one sampleset
 * \param gridScores		Output: The scores for each state on the grid
 *
 */
__global__ void
BCM_PosMeasML(const dsp::utils::state_t<double> *satStates, const cufftDoubleComplex *codeScores,
		const double *centerPt, const dsp::utils::statePosManifold_t<double> *grid, const double *enu2ecefMat,
		const double *codeFreq, const int *cpRefTOW, const int *cpElapsedEnd, const int *cpRef, const double *codePhase,
		const double *txTime, const double rxTime, const int numChan,
		const double fs, const int numSamps,
		double *gridScores)
{
	// Thread indexing variables
	unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int stride = blockDim.x * gridDim.x;


    dsp::utils::state_t<double> currPt;	// The grid point to be evaluated
    double currScore;					// The score of the current point to be evaluated


    double satPosTransmitTime;
    dsp::utils::state_t<double> currPtSatState;


    // Back-calculation values for finding scores
    double bc_los[3];
    double bc_range;
    double bc_pseudorange;
    double bc_txTime;
    double bc_codeFracDiff;
    double bc_rc;
    double bc_rangetime;
    double bc_rc0;
    double bc_rc0_idx_base;
    double bc_rc0_idx;
    int bc_rc0_cidx;
    int bc_rc0_fidx;
    cufftDoubleComplex bc_rc0_corr;


    int currTIdx;

    int test = 0;


    // Determine scores for each candidate point
    while (i < BCM_GRID_SIZE_d[0]) {

    	// Reset the score holder for this thread
    	currScore = 0;

    	// Find the ECEF position of this candidate point
    	currPt.x = enu2ecefMat[0]*grid[i].x + enu2ecefMat[1]*grid[i].y + enu2ecefMat[2]*grid[i].z + centerPt[0];
    	currPt.y = enu2ecefMat[3]*grid[i].x + enu2ecefMat[4]*grid[i].y + enu2ecefMat[5]*grid[i].z + centerPt[1];
    	currPt.z = enu2ecefMat[6]*grid[i].x + enu2ecefMat[7]*grid[i].y + enu2ecefMat[8]*grid[i].z + centerPt[2];
    	currPt.delta_t = grid[i].delta_t + centerPt[3]; // Note that this is distance (delta_t times C)

    	// Determine what time state we're looking at
		currTIdx = POSMOD(i, BCM_GRID_DIMS_d[3]);


    	// Iterate over every satellite with valid ephems
    	currScore = 0;
    	for (int currChan = 0; currChan < numChan; currChan++) {

    		// Get the pre-calculated satellite state
    		// For now, just get the satellite state for 0 delta_t grid bias
    		currPtSatState = satStates[currChan*BCM_GRID_DIMS_d[3] + (BCM_GRID_DIMS_d[3]/2)];


    		// Back-calculate the channel parameters
    		bc_los[0] = currPtSatState.x - currPt.x;
    		bc_los[1] = currPtSatState.y - currPt.y;
    		bc_los[2] = currPtSatState.z - currPt.z;
    		bc_range = norm(3, bc_los); // Might as well find range this way since the unit vector is needed.
    		bc_pseudorange = bc_range - CONST_C * currPtSatState.delta_t + currPt.delta_t;
    		bc_txTime = rxTime - bc_pseudorange/CONST_C;
    		bc_codeFracDiff = bc_txTime - cpRefTOW[currChan] - ((cpElapsedEnd[currChan] - cpRef[currChan]) * CONST_T_CA);
    		bc_rc = bc_codeFracDiff * CONST_F_CA;


			// Find the code score index
    		bc_rc0 = bc_rc - codePhase[currChan];
			bc_rc0_idx_base = (fs/codeFreq[currChan])*(-bc_rc0) + numSamps/2.0;

			// Make sure the code score index is within the range of this PRN
			// (Note: strictly less/greater so cel/flr stays within range)
			if (bc_rc0_idx_base < numSamps && bc_rc0_idx_base > 0) {
				// Also, add in the offset corresponding to this prn
				bc_rc0_idx = bc_rc0_idx_base + (numSamps*currChan);
				bc_rc0_fidx = floor(bc_rc0_idx);
				bc_rc0_cidx = floor(bc_rc0_idx+1);
			}
			else {
				//TODO: return error
				test = -4;
			}

			// Interpolate the score and add to the running total for this point
			// TODO: Do these tests need to be here? Had warp_illegal_address without them...
			cufftDoubleComplex test1 = codeScores[bc_rc0_cidx];
			cufftDoubleComplex test2 = codeScores[bc_rc0_fidx];
			cufftDoubleComplex test3 = BCM_ComplexScale(test1, bc_rc0_idx-bc_rc0_fidx);
			cufftDoubleComplex test4 = BCM_ComplexScale(test2, bc_rc0_cidx-bc_rc0_idx);
			bc_rc0_corr = cuCadd(test3, test4);
			//bc_rc0_corr = cuCadd(BCM_ComplexScale(codeScores[bc_rc0_cidx], bc_rc0_idx-bc_rc0_fidx),
			//					 BCM_ComplexScale(codeScores[bc_rc0_fidx], bc_rc0_cidx-bc_rc0_idx));

			currScore += pow(cuCabs(bc_rc0_corr), BCM_LPOWER_d);

    	} // Done iterating over all PRNs

    	// Copy score so the manifold can be measured (later)
    	gridScores[i] = currScore;

    	i += stride;
    }



}




/** \brief Compute the velocity-domain measurement as the highest score state from the manifold grid
 *
 * Note: Expecting a 1D kernel launch
 *
 * Note: SatPos computes the relativistic corrections and factors them in to the txTime provided.
 * So, when looking up a satellite position by time, it should be looked up WITHOUT application
 * of relativistic corrections (batch SatPos effectively calculates PyGNSS's
 * satellite_clock_correction immediately followed by locate_satellite).
 * However, when computing back-calculation equations, the relativistic corrections do need to be
 * applied appropriately. These can be retrieved from the delta_t and delta_t_dot elements of
 * a given satellite state.
 *
 *
 * \param satStates 		Ptr to the state of each satellite
 * \param carrScores 		Ptr to the array of carrier scores (1D array of scores, with tracked PRNs adjacent)
 * \param centerPt   		Ptr to the array of the ECEF location by which the template ENU grid should be
 *        		        	translated and rotated aka X_t|t-1 EKF estimate
 * \param grid	   			The template ENU grid
 * \param enu2ecefMat		Ptr to a rotation matrix from ENU (grid domain) to ECEF (satellite domain)
 * \param carrFreq   		Ptr to array of the carrFreqs of the tracked channels (used to find code scores idx)
 * \param txTime			Ptr to array of the transmission times from sat as calculated by the chan params
 * \param rxTime			The current time of the receiver's internal clock (no delta-t applied)
 * \param numChan			The number of satellites being tracked
 * \param fs				The sampling frequency of the sample set (used to find code scores idx)
 * \param numfftPts			The number of scores per channel
 * \param gridScores		Output: The scores for each state on the grid
 *
 */
__global__ void
BCM_VelMeasML(const dsp::utils::state_t<double> *satStates, const cufftDoubleComplex *carrScores,
		const double *centerPt, const dsp::utils::stateVelManifold_t<double> *grid, const double *enu2ecefMat,
		const double *carrFreq, const double *txTime, const double rxTime, const int numChan,
		const double fs, const int numfftPts, const int *dopplerSign,
		double *gridScores)
{
	// Thread indexing variables
	unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int stride = blockDim.x * gridDim.x;


    dsp::utils::stateVelManifold_t<double> currPt;	// The grid point to be evaluated
    double currScore;					// The score of the current point to be evaluated

    dsp::utils::state_t<double> currPtSatState;

    // Back-calculation values for finding scores
    double currPtVelECI[4];
    double bc_los[3];
    double bc_range;
    double bc_losrangerate;
    double bc_fi;
    double bc_fi0;
    double bc_fi0_idx_base;
    double bc_fi0_idx;
    int bc_fi0_fidx;
    int bc_fi0_cidx;
    cufftDoubleComplex bc_fi0_fft;

    int test = 0;
    dsp::utils::weightState_t<double> testState = {};


    // Determine scores for each candidate point
    while (i < BCM_GRID_SIZE_d[1]) {

    	// Reset the score holder for this thread
    	currScore = 0;

    	// Find the ECEF position of this candidate point (note: this is actually a velocity state)
    	currPt.x_dot = enu2ecefMat[0]*grid[i].x_dot + enu2ecefMat[1]*grid[i].y_dot + enu2ecefMat[2]*grid[i].z_dot + centerPt[4];
    	currPt.y_dot = enu2ecefMat[3]*grid[i].x_dot + enu2ecefMat[4]*grid[i].y_dot + enu2ecefMat[5]*grid[i].z_dot + centerPt[5];
    	currPt.z_dot = enu2ecefMat[6]*grid[i].x_dot + enu2ecefMat[7]*grid[i].y_dot + enu2ecefMat[8]*grid[i].z_dot + centerPt[6];
    	currPt.delta_t_dot = grid[i].delta_t_dot + centerPt[7];


    	// Iterate over every satellite with valid ephems
    	currScore = 0;
    	for (int currChan = 0; currChan < numChan; currChan++) {

    		// Assuming symmetric time grid, take the middle point in the set for this channel
    		currPtSatState = satStates[currChan*BCM_GRID_DIMS_d[3] + (BCM_GRID_DIMS_d[3]/2)];

    		// Also need to rotate the point's velocity into the inertial frame
    		currPtVelECI[0] = currPt.x_dot - CONST_OEDot*centerPt[1];
    		currPtVelECI[1] = currPt.y_dot + CONST_OEDot*centerPt[0];
    		currPtVelECI[2] = currPt.z_dot;
    		currPtVelECI[3] = currPt.delta_t_dot;

			// Find carrier frequency (using centerPt since position states are fixed in vel manifold)
			bc_los[0] = currPtSatState.x - centerPt[0];
			bc_los[1] = currPtSatState.y - centerPt[1];
			bc_los[2] = currPtSatState.z - centerPt[2];
			bc_range = norm(3, bc_los); // Might as well find range this way since the unit vector is needed.
			bc_losrangerate = 	((bc_los[0]/bc_range)*(currPtVelECI[0]-currPtSatState.x_dot)) +
								((bc_los[1]/bc_range)*(currPtVelECI[1]-currPtSatState.y_dot)) +
								((bc_los[2]/bc_range)*(currPtVelECI[2]-currPtSatState.z_dot));
			bc_fi = CONST_F_L1 * ((bc_losrangerate - currPtVelECI[3])/CONST_C + currPtSatState.delta_t_dot) / (*dopplerSign);


			bc_fi0 = bc_fi - carrFreq[currChan];
			// As opposed to bc_rc0, which would be negative here, based on the way fi is defined,
			// bc_fi0 should be positive here
			bc_fi0_idx_base = (numfftPts/fs)*(bc_fi0) + numfftPts/2.0;

			// Make sure the code score index is within the range of this PRN
			// (Note: strictly less/greater so cel/flr stays within range)
			if (bc_fi0_idx_base < numfftPts && bc_fi0_idx_base > 0) {
				// Also, add in the offset corresponding to this prn
				bc_fi0_idx = bc_fi0_idx_base + (numfftPts*currChan);
				bc_fi0_fidx = floor(bc_fi0_idx);
				bc_fi0_cidx = floor(bc_fi0_idx+1);
			}
			else {
				//TODO: return error
				test = -4;
			}

			// Interpolate the score and add to the running total for this point
			bc_fi0_fft = cuCadd(BCM_ComplexScale(carrScores[bc_fi0_cidx], bc_fi0_idx-bc_fi0_fidx),
								BCM_ComplexScale(carrScores[bc_fi0_fidx], bc_fi0_cidx-bc_fi0_idx));
			currScore += pow(cuCabs(bc_fi0_fft), BCM_LPOWER_d);

    	} // Done iterating over all channels

    	// Copy score to be accessible for display (note: this is only for display, unnecessary for functionality)
    	gridScores[i] = currScore;

    	i += stride;
    }
}



/** \brief Rotate the ENU ML state to ECEF and generate the covariance matrix
 *
 * \param stateMLIdx 		The index on the grid of the ML state
 * \param centerPt 			The state on which the grid is referenced
 * \param grid		   		The grid used in this measurement
 * \param enu2ecefMat		Ptr to a rotation matrix from ENU (grid domain) to ECEF (satellite domain)
 * \param zVal				Ptr to the array holding the measurement for this iteration
 * \param RVal				Ptr to the array holding the covariance of the measurement for this iteration
 *
 */
__global__ void
BCM_MakePosMeas(const int stateMLIdx,
		const double *centerPt, const dsp::utils::statePosManifold_t<double> *grid,
		const double *enu2ecefMat, double *zVal, double *RVal) {

	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int stride = blockDim.x * gridDim.x;

	int test = 0;

	while (i < 36) {
		switch (i) {
			case 32:
				zVal[0] = enu2ecefMat[0]*grid[stateMLIdx].x + enu2ecefMat[1]*grid[stateMLIdx].y + enu2ecefMat[2]*grid[stateMLIdx].z + centerPt[0];
				break;
			case 33:
				zVal[1] = enu2ecefMat[3]*grid[stateMLIdx].x + enu2ecefMat[4]*grid[stateMLIdx].y + enu2ecefMat[5]*grid[stateMLIdx].z + centerPt[1];
				break;
			case 34:
				zVal[2] = enu2ecefMat[6]*grid[stateMLIdx].x + enu2ecefMat[7]*grid[stateMLIdx].y + enu2ecefMat[8]*grid[stateMLIdx].z + centerPt[2];
				break;
			case 35:
				zVal[3] = grid[stateMLIdx].delta_t + centerPt[3];
				break;

			// Set the covariance matrix to identity (for now)
			case (0):
			case (9):
			case (18):
			case (27):
				RVal[i] = 1;
				break;
			default:
				RVal[i] = 0;
				break;
		}

		i += stride;
	}
}



/** \brief Rotate the ENU ML state to ECEF and generate the covariance matrix
 *
 * \param stateMLIdx 		The index on the grid of the ML state
 * \param centerPt 			The state on which the grid is referenced
 * \param grid		   		The grid used in this measurement
 * \param enu2ecefMat		Ptr to a rotation matrix from ENU (grid domain) to ECEF (satellite domain)
 * \param zVal				Ptr to the array holding the measurement for this iteration
 * \param RVal				Ptr to the array holding the covariance of the measurement for this iteration
 *
 */
__global__ void
BCM_MakeVelMeas(const int stateMLIdx,
		const double *centerPt, const dsp::utils::stateVelManifold_t<double> *grid,
		const double *enu2ecefMat,
		double *zVal, double *RVal) {

	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int stride = blockDim.x * gridDim.x;

	while (i < 36) {
		switch (i) {
			case 32:
				zVal[4] = enu2ecefMat[0]*grid[stateMLIdx].x_dot + enu2ecefMat[1]*grid[stateMLIdx].y_dot + enu2ecefMat[2]*grid[stateMLIdx].z_dot + centerPt[4];
				break;
			case 33:
				zVal[5] = enu2ecefMat[3]*grid[stateMLIdx].x_dot + enu2ecefMat[4]*grid[stateMLIdx].y_dot + enu2ecefMat[5]*grid[stateMLIdx].z_dot + centerPt[5];
				break;
			case 34:
				zVal[6] = enu2ecefMat[6]*grid[stateMLIdx].x_dot + enu2ecefMat[7]*grid[stateMLIdx].y_dot + enu2ecefMat[8]*grid[stateMLIdx].z_dot + centerPt[6];
				break;
			case 35:
				zVal[7] = grid[stateMLIdx].delta_t_dot + centerPt[7];
				break;

			// Set the covariance matrix to identity (for now)
			case (0+4):
			case (9+4):
			case (18+4):
			case (27+4):
				RVal[i+32] = 1;
				break;
			default:
				RVal[i+32] = 0;
				break;
		}

		i += stride;
	}
}





/** \brief Debug function to make viewing device memory easier. Placing a breakpoint in a kernel doesn't guarantee thread safety.
 * 		   This function makes no modifications to data, so run it on a stream to view device memory at that point.
 *
 * \param *arr		The double array you want to view
 * \param arrLen	The size of the array you want to view
 *
 */
__global__ void
BCM_valInspector(double *arr, int arrLen)
{
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int stride = blockDim.x * gridDim.x;

    while (i<arrLen) {

    	if (i == 0) {
    		double test = arr[i];
    	}

    	i += stride;
    }
}


/** \brief Debug function to make viewing device memory easier. Placing a breakpoint in a kernel doesn't guarantee thread safety.
 * 		   This function makes no modifications to data, so run it on a stream to view device memory at that point.
 *
 * \param *arr		The double array you want to view
 * \param arrLen	The size of the array you want to view
 * \param calledby	Give this a unique value to identify where you are in the stream
 *
 */
__global__ void
BCM_valInspector2(double *arr, int arrLen, int calledby)
{
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int stride = blockDim.x * gridDim.x;

    int test;

    while (i<arrLen) {

        if (i == 0) {
        	test = 1;
        }

    	i += stride;
    }
}



/** \brief Debug function to make viewing device memory easier. Placing a breakpoint in a kernel doesn't guarantee thread safety.
 * 		   This function makes no modifications to data, so run it on a stream to view device memory at that point.
 *
 * \param gridPosScores	Ptr to array of grid scores
 * \param gridVelScores	Ptr to array of grid scores
 * \param posGridMLIdx	Index of ML state on grid
 * \param velGridMLIdx	Index of ML state on grid
 * \param posGridLocs	The dx's of the grid
 * \param velGridLocs	The dx's of the grid
 * \param codePhase		Code Phase
 * \param codeFrequency	Code Frequency
 * \param carrPhase		Carrier Phase
 * \param carrFreq		Carrier Frequency
 * \param arrLen		Number of elements to iterate over
 *
 */
__global__ void
BCM_valInspector3(double *gridPosScores, double *gridVelScores, int posGridMLIdx, int velGridMLIdx, dsp::utils::statePosManifold_t<double> *posGridLocs, dsp::utils::stateVelManifold_t<double> *velGridLocs, double *codePhase, double *codeFrequency, double *carrPhase, double *carrFreq, int arrLen)
{
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int stride = blockDim.x * gridDim.x;

    //int currChan;
    //int currSamp;

    int test1;
    int test2;
    int test3;
    int test4;
    int test5;
    int test6;
    int test7;
    int test8;
    int test9;
    int test10;


    while (i<arrLen) {

        if (i == 0) {
        	test1 = codePhase[i];
        	test2 = codeFrequency[i];
        	test3 = carrPhase[i];
        	test4 = carrFreq[i];
        	test5 = gridPosScores[i];
        	test6 = gridVelScores[i];
        	test7 = posGridMLIdx;
        	test8 = velGridMLIdx;
        	test9 = posGridLocs[i].x;
        	test10 = velGridLocs[i].x_dot;
        	// Break point here
        }

    	i += stride;
    }
}


/** \brief Debug function to make viewing device memory easier. Placing a breakpoint in a kernel doesn't guarantee thread safety.
 * 		   This function makes no modifications to data, so run it on a stream to view device memory at that point.
 *
 * \param *posGridLocs	The statePosManifold array you want to view
 * \param arrLen		The size of the array you want to view
 *
 */
__global__ void
BCM_valInspector4(dsp::utils::statePosManifold_t<double> *posGridLocs, int arrLen)
{
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int stride = blockDim.x * gridDim.x;

    double test = 0;


    while (i<arrLen) {

    	if (i == 0) {
    		test = posGridLocs[i].x;
    	}

    	i += stride;
    }

}


/** \brief Debug function to make viewing device memory easier. Placing a breakpoint in a kernel doesn't guarantee thread safety.
 * 		   This function makes no modifications to data, so run it on a stream to view device memory at that point.
 *
 * \param arr1		The array to subtract elementwise from arr2
 * \param arr2 		The reference array which arr1 is compared to
 * \param arrSize	The number of elements in each array
 *
 */
__global__ void
BCM_valInspectorDiff(double *arr1, double *arr2, int arrSize) {
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

	double diff[8];

	int test = 0;

	if (i == 0) {
		for (int idx = 0; idx < arrSize; idx++) {
			diff[idx] = arr2[idx] - arr1[idx];
		}

		test = 1;
	}
}





/*
 *
 * HOST FUNCTIONS
 *
 */

dsp::BatchCorrManifold::BatchCorrManifold() {

    ModuleName = "BatchCorrManifold";
    AllocateInputs(19);
    AllocateOutputs(4);

    Started = 0;

    /**
     * INPUT CONFIGURATION
     */
    // Configure inputs
    // Technically, the defines in SatPos can be used instead of SatPosArrSize,
    // but it is passed as an input for future accessibility purposes
    ConfigExpectedInput(0, "CodeScores", 		UNDEFINED_t, 	VALUE_CMPX, 	VECTORLENGTH_ANY);
    ConfigExpectedInput(1, "CarrScores", 		UNDEFINED_t, 	VALUE_CMPX, 	VECTORLENGTH_ANY);
    ConfigExpectedInput(2, "xCurrkk1", 			DOUBLE_t, 		STATE, 			VECTORLENGTH_ANY);
    ConfigExpectedInput(3, "txTime",			DOUBLE_t,		VALUE,			VECTORLENGTH_ANY);
    ConfigExpectedInput(4, "SatStates",			DOUBLE_t,		STATE,			VECTORLENGTH_ANY);
    ConfigExpectedInput(5, "rxTime", 			DOUBLE_t, 		VALUE, 			1);
    ConfigExpectedInput(6, "SampleLength", 		DOUBLE_t, 		VALUE, 			1);
    ConfigExpectedInput(7, "SamplingFrequency", DOUBLE_t, 		FREQUENCY_HZ, 	1);
    ConfigExpectedInput(8, "CodeFrequency",		DOUBLE_t,		FREQUENCY_HZ,	VECTORLENGTH_ANY);
    ConfigExpectedInput(9, "CarrierFrequency",	DOUBLE_t,		FREQUENCY_HZ,	VECTORLENGTH_ANY);
    ConfigExpectedInput(10, "DopplerSign",		INT_t,			VALUE,			1);
    ConfigExpectedInput(11, "NumFFTPoints",		INT_t,			VALUE,			1);
    ConfigExpectedInput(12, "ENU2ECEFMat",		DOUBLE_t,		VALUE,			9);
    ConfigExpectedInput(13, "SatStatesOld",		DOUBLE_t,		STATE,			VECTORLENGTH_ANY);
    ConfigExpectedInput(14, "CodePhase",		DOUBLE_t,		VALUE,			VECTORLENGTH_ANY);
    ConfigExpectedInput(15, "CarrierPhase",		DOUBLE_t,		VALUE,			VECTORLENGTH_ANY);
    ConfigExpectedInput(16, "cpRefTOW",			INT_t,			VALUE,			VECTORLENGTH_ANY);
    ConfigExpectedInput(17, "cpElapsedEnd",		INT_t,			VALUE,			VECTORLENGTH_ANY);
    ConfigExpectedInput(18, "cpRef",			INT_t,			VALUE,			VECTORLENGTH_ANY);


    // TODO: Consider changing GridDimSize and GridDimSpacing to arrays for custom grid sizing
    //       (would probably need to parse a string to ints in DPInit or somewhere)
    // Get grid size
    InsertParam("PosGridDimSize", 		(void*)&posGridDimSizeParam, 		INT_t, 		sizeof(int), 	sizeof(int));
    InsertParam("VelGridDimSize", 		(void*)&velGridDimSizeParam, 		INT_t, 		sizeof(int), 	sizeof(int));
    InsertParam("GridDimSpacing", 	(void*)&gridDimSpacingParam, 	FLOAT_t, 	sizeof(float), 	sizeof(float));
    InsertParam("GridType", 		(void*)&gridTypeParam,			INT_t,
    		sizeof(dsp::utils::ManifoldGridTypes), sizeof(dsp::utils::ManifoldGridTypes));
    InsertParam("LPower", 			(void*)&LPower, 				INT_t, 		sizeof(int), 	sizeof(int));
    InsertParam("GridLogFileName", (void*)&Filename, CHAR_t, FilenameCapacity, 0);
    InsertParam("LoadPosGrid",		(void*)&loadPosGrid, BOOL_t, sizeof(bool), sizeof(bool));
    InsertParam("LoadPosGridFilename", (void*)&loadPosGridFilename, CHAR_t, FilenameCapacity, 0);


    // Configure outputs
    ConfigOutput(0, "zVal", 	DOUBLE_t, STATE, 		CUDA_DEVICE, VECTORLENGTH_ANY, NULL, 0);
    ConfigOutput(1, "RVal", 	DOUBLE_t, COVARIANCE, 	CUDA_DEVICE, VECTORLENGTH_ANY, NULL, 0);
    ConfigOutput(2, "TimeGrid", DOUBLE_t, VALUE,		CUDA_DEVICE, VECTORLENGTH_ANY, NULL, 0);
    ConfigOutput(3, "PosScores", DOUBLE_t, GRID,		CUDA_DEVICE, VECTORLENGTH_ANY, NULL, 0);

    std::clog << "[" << ModuleName << "] Constructed" << std::endl;
}


dsp::BatchCorrManifold::~BatchCorrManifold() {
	if (Started) Stop();
    delete [] inputs;
    delete [] outputs;
    delete [] expectedInputs;
}


int
dsp::BatchCorrManifold::Start(void* cuFlowStream) {

	// Check module status and report accordingly
	if (Started) {
        std::clog << "[" << ModuleName << "] Start: Already Started." << std::endl;
        return 0;
    }
    std::clog << "[" << ModuleName << "] Starting ... " << std::flush;


    // TODO: (to be updated accordingly once more sophisticated grid parameters can be provided)
    // Using the inputs, set up to store grid parameters
    // If the grid is even-sized in a given dimension, the center point will be the smaller of the two "middle points" wrt the buffer size
    gridDims[0] = gridDims[1] = gridDims[2] = gridDims[3] = posGridDimSizeParam;
    gridDims[4] = gridDims[5] = gridDims[6] = gridDims[7] = velGridDimSizeParam;
    for (int i = 0; i < 2*BCM_NUM_GRID_STATES; i++) {
    	gridHalfIdx[i] = (gridDims[i] - 1)/2;
    	gridDimSpacingArr[i] = gridDimSpacingParam;
    }

    // *** HACK TO LOCALLY CONFIGURE THE GRID HERE ***
    // TODO: Remove once more sophisticated grid parameters setup is done
    // Example configuration is commented below ------->
    //gridDimSpacingArr[0] = gridDimSpacingArr[1] = 1.5;
    //gridDimSpacingArr[2] = 3.0;
    //gridDimSpacingArr[3] = 12.5;
    //gridDimSpacingArr[4] = gridDimSpacingArr[5] = gridDimSpacingArr[0] / 10.0;
    //gridDimSpacingArr[6] = gridDimSpacingArr[2] / 10.0;
    //gridDimSpacingArr[7] = 0.25;


	// Set the CUDA Stream for the GPU operations
    cuStream = (cudaStream_t*)cuFlowStream;


    // Check grid size
    gridSize[0] = gridDims[0] * gridDims[1] * gridDims[2] * gridDims[3];
    gridSize[1] = gridDims[4] * gridDims[5] * gridDims[6] * gridDims[7];
    if(gridSize[0]+gridSize[1] > BCM_MAX_GRID_SIZE) {
    	std::clog << "[" << ModuleName << "] Manifolds specified exceed maximum size of " << BCM_MAX_GRID_SIZE <<  " points."
    			<< std::endl << "    Position manifold: " << gridSize[0] << " points."
    			<< std::endl << "    Velocity manifold: " << gridSize[1] << " points."
    			<< std::endl << "This size restriction is based on expected GPU memory usage and the limit of the Jetson TX2."
    			<< std::endl << "Either reduce the grid sizes so their sum is below BCM_MAX_GRID_SIZE (batchcorrmanifold.h)"
    			<< std::endl << "or make sure you have enough memory before updating BCM_MAX_GRID_SIZE!" << std::endl;
    	return -1;
    }

    gridPosScores_h = new double[gridSize[0]];

    // Set up streams for processing
	cuCheck(cudaStreamCreate(&posStream));
	cuCheck(cudaStreamCreate(&velStream));


    // Set device constants so these values don't have to be calculated/computed each iteration
    cuCheckMSt(cudaMemcpyToSymbol(BCM_GRID_DIMS_d, 			&gridDims, 					sizeof(int)*2*BCM_NUM_GRID_STATES, 	0, cudaMemcpyHostToDevice));
    cuCheckMSt(cudaMemcpyToSymbol(BCM_HALF_IDX_d, 			&gridHalfIdx, 				sizeof(int)*2*BCM_NUM_GRID_STATES, 	0, cudaMemcpyHostToDevice));
    cuCheckMSt(cudaMemcpyToSymbol(BCM_GRID_SIZE_d, 			&gridSize, 					sizeof(int)*2, 						0, cudaMemcpyHostToDevice));
    cuCheckMSt(cudaMemcpyToSymbol(BCM_LPOWER_d, 			&LPower, 					sizeof(int), 						0, cudaMemcpyHostToDevice));


    // Allocate space
    // Create position grid/scores (Indexing is grid[xidx][yidx][zidx][tidx])
    cuCheckMSt(cudaMalloc((void**)&gridPosScores_d, sizeof(double) * gridSize[0]));
    cuCheckMSt(cudaMalloc((void**)&gridPosLocs_d, 	sizeof(dsp::utils::statePosManifold_t<double>) * gridSize[0]));
    // Create velocity grid/scores (Indexing is grid[xdotidx][ydotidx][zdotidx][tdotidx])
    cuCheckMSt(cudaMalloc((void**)&gridVelScores_d, sizeof(double) * gridSize[1]));
    cuCheckMSt(cudaMalloc((void**)&gridVelLocs_d, 	sizeof(dsp::utils::stateVelManifold_t<double>) * gridSize[1]));
    // Create output states
    cuCheckMSt(cudaMalloc((void**)&zVal_d, sizeof(dsp::utils::state_t<double>)));
    cuCheckMSt(cudaMalloc((void**)&RVal_d, sizeof(double) * (2*BCM_NUM_GRID_STATES) * (2*BCM_NUM_GRID_STATES)));
    // Store grid properties on device
    cuCheckMSt(cudaMalloc((void**)&gridDimSpacing_d, sizeof(double)*(2*BCM_NUM_GRID_STATES)));
    cuCheckMSt(cudaMemcpyAsync(gridDimSpacing_d, gridDimSpacingArr, sizeof(double)*(2*BCM_NUM_GRID_STATES), cudaMemcpyHostToDevice, *cuStream));
    // Allocate an (effectively) 2D array for the reduction sum
    cuCheckMSt(cudaMalloc((void**)&weightedPosStates_d, sizeof(dsp::utils::weightState_t<double>)*threadsPerMeasRedBlock));
    cuCheckMSt(cudaMalloc((void**)&weightedVelStates_d, sizeof(dsp::utils::weightState_t<double>)*threadsPerMeasRedBlock));
    // Allocate the time grid for pre-processing
    cuCheckMSt(cudaMalloc((void**)&timeGrid_d, sizeof(double)*gridDims[3]));


    // Now that the space is allocated, the output pointers can be assigned
    outputs[0].Data = zVal_d;
    outputs[0].VectorLength = 2*BCM_NUM_GRID_STATES;
    outputs[1].Data = RVal_d;
    outputs[1].VectorLength = (2*BCM_NUM_GRID_STATES) * (2*BCM_NUM_GRID_STATES);
    outputs[2].Data = timeGrid_d;
    outputs[2].VectorLength = gridDims[3];
    outputs[3].Data = (void*)gridPosScores_d;
    outputs[3].VectorLength = (unsigned short)gridSize[0];

    // Type cast and store the input pointers for accessibility
    codeScores_d	= ((cufftDoubleComplex*)(inputs[0]->Data));
    carrScores_d	= ((cufftDoubleComplex*)(inputs[1]->Data));
    // Inputs 3 and 4 get copied to device constants; don't need pointers to them

    // Get the inputs accessible
    numfftPointsPtr	= ((int*)(inputs[11]->Data));


    // And generate the manifolds
    // (Hack -- default initialize the grid so timeGrid_d gets created, then overwrite it if we're loading from file.
    //  This hack is acceptable because the only kind of grid we want to load from file is the RNGrid,
    //  which is random in time offsets, making timeGrid_d not useful. If a structured grid is to be loaded from file,
    //  some better way to handle timeGrid_d should be considered.)
    BCM_InitPosGrid<<<1, 512, 0, *cuStream>>>(gridPosLocs_d, gridTypeParam, gridDimSpacing_d, timeGrid_d);
    if (loadPosGrid) {
    	dsp::utils::statePosManifold_t<double> * loadGridVals = new dsp::utils::statePosManifold_t<double>[gridSize[0]];
    	const char* tok;
    	int curGridIdx = 0;

        FILE* loadGridFile = fopen(loadPosGridFilename, "r");
        if (loadGridFile == NULL) {
        	std::clog << "[" << ModuleName << "]" << " Open loadGridFile failed: " << loadPosGridFilename << std::endl;
        	return -1;
        }

        char curLine[1024];
        while (fgets(curLine, 1024, loadGridFile)) {
            const char* val1 = strtok(curLine, ",");
			loadGridVals[curGridIdx].x = atof(val1);
			const char* val2 = strtok(NULL, ",");
			loadGridVals[curGridIdx].y = atof(val2);
			const char* val3 = strtok(NULL, ",");
			loadGridVals[curGridIdx].z = atof(val3);
			const char* val4 = strtok(NULL, "\r");
			loadGridVals[curGridIdx].delta_t = atof(val4);
			curGridIdx++;
        }

        // Copy the loaded grid offsets to device memory
        cuCheckMSt(cudaMemcpyAsync(gridPosLocs_d, loadGridVals, sizeof(dsp::utils::statePosManifold_t<double>)*(gridSize[0]), cudaMemcpyHostToDevice, *cuStream));
    }
    BCM_InitVelGrid<<<1, 512, 0, *cuStream>>>(gridVelLocs_d, gridTypeParam, gridDimSpacing_d);


    BCM_valInspector4<<<1,64,0, *cuStream>>>(gridPosLocs_d, 1);

    // Make sure all GPU tasks have completed before continuing
    cuCheckMSt(cudaStreamSynchronize(*cuStream));


    // Signifies that the next call to update() will be the first after start()
    Started = 1;

    std::clog << "Started." << std::endl;
    return 0;
}


int
dsp::BatchCorrManifold::Stop() {

    int ret = 0;
    if (Started == 0) {
        std::clog << "[" << ModuleName << "] Stop: Wasn't Started." << std::endl;
        return 0;
    }
    std::clog << "[" << ModuleName << "] Stopping ... " << std::flush;

    // Destroy ancillary stream
    cuCheckV(cudaStreamDestroy(posStream));
    cuCheckV(cudaStreamDestroy(velStream));


    // Free device memory
    cuCheckMSp(cudaFree((void*)gridPosScores_d));
    cuCheckMSp(cudaFree((void*)gridVelScores_d));
    cuCheckMSp(cudaFree((void*)gridPosLocs_d));
    cuCheckMSp(cudaFree((void*)gridVelLocs_d));
    cuCheckMSp(cudaFree((void*)zVal_d));
    cuCheckMSp(cudaFree((void*)RVal_d));
    cuCheckMSp(cudaFree((void*)gridDimSpacing_d));
    cuCheckMSp(cudaFree((void*)timeGrid_d));
    cuCheckMSp(cudaFree((void*)weightedPosStates_d));
    cuCheckMSp(cudaFree((void*)weightedVelStates_d));


    Started = 0;
    std::clog << "Stopped." << std::endl;

    return ret;
}


int
dsp::BatchCorrManifold::Update(void* cuFlowStream) {

    if (Started == 0) {
        std::cerr << "[" << ModuleName
                  << "] Error: Update() Failed due to SatPos not initialized"
                  << std::endl;
        return -1;
    }


    if (FirstUpdate) {

    	// Get the inputs from the downstream modules
    	numChan 		= inputs[3]->VectorLength;
        txTimePtr_d 	= ((double*)(inputs[3]->Data));
        rxTimePtr 		= ((double*)(inputs[5]->Data));
        sampleLengthPtr = ((double*)(inputs[6]->Data));
        samplingFrequencyPtr = ((double*)(inputs[7]->Data));
        codeFrequency_d = ((double*)(inputs[8]->Data));
        carrFrequency_d = ((double*)(inputs[9]->Data));
        dopplerSign_d	= ((int*)(inputs[10]->Data));
        satStates_d 	= ((dsp::utils::state_t<double>*)(inputs[4]->Data));
        enu2ecefMat_d 	= ((double*)(inputs[12]->Data));
        satStatesOld_d  = ((dsp::utils::state_t<double>*)(inputs[13]->Data));
        codePhase_d = ((double*)(inputs[14]->Data));
        carrPhase_d = ((double*)(inputs[15]->Data));
        cpRefTOW_d = ((int*)(inputs[16]->Data));
        cpElapsedEnd_d = ((int*)(inputs[17]->Data));
        cpRef_d = ((int*)(inputs[18]->Data));

        FirstUpdate = 0;
    }

    // Update the rxTime to consider this set of samples
    numSamps = (*sampleLengthPtr) * (*samplingFrequencyPtr);
    currSamplingFreq = *samplingFrequencyPtr;

    // Get the center point for this iteration's manifold
    xCurr_d = (double*)inputs[2]->Data;



    // Evaluate the manifolds

    // Method 1: Weighted-average method
    /*
    BCM_PosMeasReduction<<<threadsPerMeasRedBlock, threadsPerPosManiBlock, sizeof(dsp::utils::weightState_t<double>)*threadsPerPosManiBlock, posStream>>>
    		(satStates_d, codeScores_d,
    		 xCurr_d, gridPosLocs_d, enu2ecefMat_d,
    		 codeFrequency_d, txTimePtr_d, *rxTimePtr, numChan,
    		 currSamplingFreq, numSamps,
    		 weightedPosStates_d);

    BCM_VelMeasReduction<<<threadsPerMeasRedBlock, threadsPerVelManiBlock, sizeof(dsp::utils::weightState_t<double>)*threadsPerVelManiBlock, velStream>>>
    		(satStates_d, carrScores_d,
    		 xCurr_d, gridVelLocs_d, enu2ecefMat_d,
    		 carrFrequency_d, txTimePtr_d, *rxTimePtr, numChan,
    		 currSamplingFreq, *numfftPointsPtr, dopplerSign_d,
    		 weightedVelStates_d);

    BCM_ReduceAndPosMeas<<<1, threadsPerMeasRedBlock, sizeof(dsp::utils::weightState_t<double>)*threadsPerMeasRedBlock , posStream>>>
    		(weightedPosStates_d, threadsPerMeasRedBlock, zVal_d, RVal_d);

    BCM_ReduceAndVelMeas<<<1, threadsPerMeasRedBlock, sizeof(dsp::utils::weightState_t<double>)*threadsPerMeasRedBlock , velStream>>>
    		(weightedVelStates_d, threadsPerMeasRedBlock, zVal_d, RVal_d);
     */


    // Method 2: Maximum likelihood (highest score state) method
    // NOTE: If you change threadsPer___ManiBlock, you need to change the shared memory allocation in the kernel itself to match!
    //       Currently hacking a compile time-known value in the kernel so weird extern allocation stuff can be ignored.
    BCM_PosMeasML<<<threadsPerMeasRedBlock, threadsPerPosManiBlock, 0, posStream>>>
    		(satStates_d, codeScores_d,
    		 xCurr_d, gridPosLocs_d, enu2ecefMat_d,
    		 codeFrequency_d, cpRefTOW_d, cpElapsedEnd_d, cpRef_d, codePhase_d,
    		 txTimePtr_d, *rxTimePtr, numChan,
    		 currSamplingFreq, numSamps,
    		 gridPosScores_d);

    BCM_VelMeasML<<<threadsPerMeasRedBlock, threadsPerVelManiBlock, 0, velStream>>>
    		(satStates_d, carrScores_d,
    		 xCurr_d, gridVelLocs_d, enu2ecefMat_d,
    		 carrFrequency_d, txTimePtr_d, *rxTimePtr, numChan,
    		 currSamplingFreq, *numfftPointsPtr, dopplerSign_d,
    		 gridVelScores_d);

    // Find the pointer idx for the manifold's largest score
    gridPosMaxPtr_d = thrust::max_element(thrust::cuda::par.on(posStream), gridPosScores_d, gridPosScores_d + gridSize[0]);
    gridVelMaxPtr_d = thrust::max_element(thrust::cuda::par.on(velStream), gridVelScores_d, gridVelScores_d + gridSize[1]);

    //BCM_valInspector<<<1,64,0, velStream>>>(gridVelScores_d, 101); // Debug

    // Process this largest score as a measurement
    BCM_MakePosMeas<<<1, 64, 0, posStream>>>((int)(gridPosMaxPtr_d - gridPosScores_d), xCurr_d, gridPosLocs_d, enu2ecefMat_d, zVal_d, RVal_d);
    BCM_MakeVelMeas<<<1, 64, 0, velStream>>>((int)(gridVelMaxPtr_d - gridVelScores_d), xCurr_d, gridVelLocs_d, enu2ecefMat_d, zVal_d, RVal_d);

    //BCM_valInspector<<<1,64,0, posStream>>>(carrFrequency_d, 102); // Debug
    //BCM_valInspector<<<1,64,0, posStream>>>(codePhase_d, 103); // Debug





    // Block on host until both manifolds have completed
    cuCheckMSt(cudaStreamSynchronize(posStream));
    cuCheckMSt(cudaStreamSynchronize(velStream));




    // Grid logging
    // USE THIS WITH CAUTION! Can become very space and time-expensive
    /*
    cuCheckM(cudaMemcpyAsync(gridPosScores_h, gridPosScores_d, sizeof(double)*gridSize[0],
             cudaMemcpyDeviceToHost, *cuStream));

	std::ofstream gridFile;
	gridFile.open(Filename);

	for (int idx = 0; idx < gridSize[0]; idx++) {
		gridFile << gridPosScores_h[idx];
		if (idx < gridSize[0] - 1) {
			gridFile << ", ";
		}
	}
	gridFile << std::endl;
	gridFile.close();
     */

    //BCM_valInspectorDiff<<<1, 64, 0, *cuStream>>>(zVal_d, xCurr_d, 8); // Debug
    cuCheckMSt(cudaStreamSynchronize(*cuStream));

	return 0;
}

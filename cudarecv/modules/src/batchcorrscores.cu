#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <cuComplex.h>
#include <cufft.h>
#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <stdint.h>
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
#include "batchcorrscores.h"
#include "errorhandler.h"
#include "auxil.h"
#include "ephhelper.h"
//#define ACQ_TEST

#define REDUCE_MAX_INITIAL_BLKSIZE  512
#define PREPROCESS_CA_FFT 0
//__constant__ cufftDoubleComplex deviceCACode[PRN_MAX][1023];

__device__ __constant__ int BCS_LCA_d;
__device__ __constant__ double BCS_FS_d;
__device__ __constant__ int BCS_S_d;
__device__ __constant__ double BCS_T_d;
__device__ __constant__ int BCS_SATPOS_ARR_SIZE_d;
__device__ __constant__ int BCS_SATPOS_NUM_ARRS_d;
__device__ __constant__ double BCS_SATPOS_BATCH_TIME_d;

/**
 * ==============================================================
 * || MODULE SUMMARY                                           ||
 * ==============================================================
 * 
 * Generates the correlation scores for the position-time and velocity-drift domains.
 * These scores are indexed by their channel parameters: code phase and carrier frequency.
 * 
 * Note: Channel and channel parameters in this context referes to the information
 * needed to follow one specific satellite:
 *      - Code phase
 *      - Carrier phase
 *      - Code frequency
 *      - Carrier frequency
 *      - Code periods elapsed since start 
 *      - Reference code period of last Z-count
 * Refer to Chapters 2 and 3 of Matthew Peretic's MS thesis
 * (http://hdl.handle.net/2142/10523) for mathematical definitions of these terms
 * 
 * Though these parameters are not being used to estimate intermediate measurements 
 * (like pseudoranges), knowing them wrt the last timestep's navigation solution 
 * allows for "batch correlation" -- specifically, calculating correlation values
 * for replicas with parameters close to the parameters of the last solution.
 * To solve for a navigation solution, then, the individual grid points just need
 * to calculate these channel parameters and interpolate between the batch-calculated
 * values.
 * 
 * Alternately, each grid point could generate its own signal replica from scratch
 * and run a reduction-sum for the length of the sample set to perform correlation.
 * This would have to be run for every point on an 8-D grid (velocity-drift must
 * be accounted for in addition to position-time to account for Doppler effects),
 * which would (likely) be both impractically time and memory-consuming. Also,
 * this would have to be calculated at a critical section of code, whereas the 
 * Batch Correlation method may be run in parallel with the satellite position calculator.
 * 
 * 
 */

/** \brief Perform addition of cufftDoubleComplex values
 *
 * Performs addition of cufftDoubleComplex values for use in thrust reductions.
 *
 * \param lhs      One cufftDoubleComplex value.
 * \param rhs      The other cufftDoubleComplex value
 *
 */
struct addCufftComplex {
	__device__ cufftDoubleComplex operator()(const cufftDoubleComplex lhs,
			const cufftDoubleComplex rhs) {
		cufftDoubleComplex temp = lhs;

		temp.x += rhs.x;
		temp.y += rhs.y;

		return temp;
	}
};



/** \brief Transposes complex values.
 *
 * \param vec Input and output complex value. Values are transposed and
 *            written back to same memory location.
 * \param len Number of complex values to transpose.
 */
__global__ void BCS_ComplexVectorTranspose(cufftDoubleComplex *vec,
		unsigned int len) {

	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int stride = blockDim.x * gridDim.x;

	while (i < len) {
		vec[i] = cuConj(vec[i]);
		i += stride;
	}
}




/** \brief Generates CA code chips
 *
 * \param prn  PRN of CA code to generate
 * \param code Buffer where code is written to.
 * \return 0 on success. -1 otherwise.
 */
__global__ void BCS_GenCACode(cufftDoubleComplex *code, int numChan) {

	// Determine current PRN from thread index
	// (Note: +1 so these values start at 1, per the convention of real PRNs.
	//  Though this +1 gets subtracted off during use, it makes it clearer to
	//  the reader what PRN sequence is being generated when inspecting a thread.
	//  Since this is only run in Start(), the added computation does not affect
	//  speed within the lööp.)
	unsigned int prn = (blockIdx.x * blockDim.x + threadIdx.x) + 1;

	if (prn < numChan) {
		// The taps to be used for generating the G2 output.  From Kaplan and
		// Hegarty.
		int8_t tap1[] = { 2, 3, 4, 5, 1, 2, 1, 2, 3, 2, 3, 5, 6, 7, 8, 9, 1, 2,
				3, 4, 5, 6, 1, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 4, 1, 2, 4 };
		int8_t tap2[] = { 6, 7, 8, 9, 9, 10, 8, 9, 10, 3, 4, 6, 7, 8, 9, 10, 4,
				5, 6, 7, 8, 9, 3, 6, 7, 8, 9, 10, 6, 7, 8, 9, 10, 10, 7, 8, 10 };

		// We will generate the CA code in terms of +/- 1.  A binary 0 is
		// represented by a +1 while a binary 1 is represented by a -1.  The initial
		// state of both registers is a binary 1 (or a -1 in our notation);
		int8_t reg1[] = { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
		int8_t reg2[] = { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };

		// Create the variables that will hold the g1 and g2 codes
		int8_t g1[1023];
		int8_t g2[1023];
		int8_t feedback1, feedback2;
		for (int i = 0; i < 1023; i++) {
			g1[i] = 0;
			g2[i] = 0;
		}

		// Perform the shifting
		for (int i = 0; i < 1023; i++) {	// 1:1023
			// Read out the correct value from the two registers
			g1[i] = reg1[9];
			g2[i] = reg2[tap1[prn - 1] - 1] * reg2[tap2[prn - 1] - 1];

			// Calculate the value to be fed back to the beginning of each shift
			// register
			feedback1 = reg1[2] * reg1[9];
			feedback2 = reg2[1] * reg2[2] * reg2[5] * reg2[7] * reg2[8]
					* reg2[9];

			// Perform the shift
			for (int j = 9; j > 0; j--) {
				reg1[j] = reg1[j - 1];
				reg2[j] = reg2[j - 1];
			}
			reg1[0] = feedback1;
			reg2[0] = feedback2;
		}

		// Multiply the two together to form the CA code
		for (int i = 0; i < 1023; i++) {
			code[(prn - 1) * 1023 + i] = make_cuDoubleComplex(
					(double) (-g1[i] * g2[i]), 0.0);
		}
	}
}

/** \brief Compute the time offset of each index
 *
 * The time index of each sample is calculated wrt the first sample
 *
 * \param timeIdcs     The buffer to place the idcs in.
 */
__global__ void BCS_GenTimeIdcs(double *timeIdcs) {

	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int stride = blockDim.x * gridDim.x;

	while (i < BCS_S_d) {
		timeIdcs[i] = (double) i / BCS_FS_d;
		// Round to the nearest 1ns to avoid garbage precision
		timeIdcs[i] = round(timeIdcs[i] * 1.0e9) / 1.0e9;
		i += stride;
	}
}

/** \brief Load CUDA kernel
 *
 * Loads raw input samples and converts it to cuFloatComplex datatype.
 *
 * \param rawFixed     Raw input samples. Expected format is interleaved
 *                     int16_t real and int16_t imaginary values.
 * \param compSamps    Buffer to store loaded and mixed samples.
 *                     Format is equivalent to bins[bin][sample]. This buffer is
 *                     expected to have capacity of bins[FreqBinsPerIter][S]
 * \param S            Number of samples per 1ms period.
 */
__global__ void BCS_Load(int16_t *rawFixed, cufftDoubleComplex *compSamps,
		const int S) {
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int stride = blockDim.x * gridDim.x;

	while (i < S) {
		// load raw sample and convert to double complex
		compSamps[i] = make_cuDoubleComplex((double) rawFixed[i * 2],
				(double) rawFixed[i * 2 + 1]);
		// stride through the operation
		i += stride;
	}
}

/** \brief Find the sample index of the next nav bit
 *
 * Find where a potential bit flip could happen with respect to the
 * start of the sample set.
 *
 * \param cpElapsed    Number of code periods since start of processing
 *                     including intialization, if applicable)
 * \param cpReference  Number of code periods since ref Z-count
 * \param codePhase    Pointer to array of code phases for each channel
 * \param codeFreq     Pointer to array of code frequencies for each channel
 * \param numChan      Number of active channels this iteration
 * \param idxNext      Pointer to array of ints of indices of next nav bit
 *
 */
__global__ void BCS_NavBitBoundary(int *cpElapsed, int *cpReference,
		double *codePhase, double *codeFreq, int numChan, int *idxNext) {
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int stride = blockDim.x * gridDim.x;

	// Since kernels are launched in multiples of 32 threads for efficiency,
	// only compute if within the range to estimate
	if (i < numChan) {
		// One nav bit is transmitted using 20 code periods,
		// so count the code periods until the next nav bit
		int cpSincePrevNavBit = POSMOD((cpElapsed[i] - cpReference[i]), 20);
		int cpToNextNavBit = 20 - cpSincePrevNavBit;

		// Scale the number of code periods to next nav bit by the
		// ratio between the sampling frequency and the estimated code frequency
		// to find the sample where the transition happens.
		idxNext[i] = std::floor((BCS_LCA_d * cpToNextNavBit - codePhase[i]) * (BCS_FS_d / codeFreq[i])) + 1;

		i += stride;
	}

}






/** \brief Generate the Doppler wipeoff signal
 *
 * Calculates the Doppler wipeoff value for each sample
 *
 * \param carrFreq	Pointer to array of carrier frequencies for each channel
 * \param carrPhase	Pointer to array of carrier phases for each channel
 * \param S			The number of samples in this iteration's sampleset
 * \param numChan	The number of channels (satellites) being tracked
 * \param timeIdc	Pointer to an array of the time difference between each sample and the first sample
 * \param doppWipe	Output array to hold the generated Doppler wipeoff
 *
 */
__global__ void BCS_ComputeDopplerWipeoff(double *carrFreq, double *carrPhase,
		int S, int numChan, double *timeIdc, cufftDoubleComplex *doppWipe) {
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int stride = blockDim.x * gridDim.x;

	int currChan; // The channel this thread is processing
	int currSamp; // The index of the sample in the sequence

	int arrLen = S * numChan;
	while (i < arrLen) {

		// Update index trackers
		currChan = i / S;
		currSamp = POSMOD(i, S);

		// Compute exp(j 2pi (freq t + phi))
		cufftDoubleComplex doppSamp;
		sincos(2 * CONST_PI	* (carrFreq[currChan] * timeIdc[currSamp] + carrPhase[currChan]), &doppSamp.y, &doppSamp.x);
		// To get exp(-j 2pi (freq t + phi)),
		// - Euler's: exp(ix) = cos(x) + jsin(x)
		// - cos is an even function, so cos(-x) = cos(x)
		// - sin is an odd function, so sin(-x) = -sin(x)
		// Thus, just need to take the complex conj
		doppWipe[i] = cuConj(doppSamp);

		// stride through the operation
		i += stride;
	}
}

/** \brief Generate the replica signal
 *
 * Calculates the replica value for each sample
 *
 * \param codeFreq		Pointer to the array of code frequencies for each channel
 * \param codePhase		Pointer to the array of code phases for each channel
 * \param S				The number of samples in this iteration's sampleset
 * \param numChan		The number of channels (satellites) being tracked
 * \param PRNs			PRN number for each tracked channel
 * \param timeIdc		Pointer to an array of the time difference between each sample and the first sample
 * \param idxNextNav	The sample index (for each channel) where the next nav bit occurs
 * \param codeChips		The full CA code sequence for each satellite, concatinated
 * \param repNoFlip		Output: the replica signal for each channel assuming that the next nav bit is the same polarity as the previous
 * \param repFlipped	Output: the replica signal for each channel assuming that the next nav bit is the opposite polarity as the previous
 *
 */
__global__ void BCS_ComputeCodeReplica(double *codeFreq, double *codePhase,
		int S, int numChan, uint8_t *PRNs, double *timeIdc, int *idxNextNav,
		cufftDoubleComplex *codeChips, cufftDoubleComplex *repNoFlip,
		cufftDoubleComplex *repFlipped) {
		//, double *repNoFlipDouble, double *repFlippedDouble) {
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int stride = blockDim.x * gridDim.x;

	int currChan; // The channel this thread is processing
	int currSamp; // The index of the sample in the sequence

	int arrLen = S * numChan;

	__syncthreads();

	while (i < arrLen) {

		// Update index trackers
		currChan = i / S;
		currSamp = POSMOD(i, S);

		// Find the index of the chip to be taken for this sample
		// - advance through the code to the corresponding prn
		// - find the chip for this prn according to where in the sequence this sample lies
		int sampIdx = (int) ((PRNs[currChan] - 1) * BCS_LCA_d)
				+ POSMOD((int)std::floor(timeIdc[currSamp]* codeFreq[currChan] + codePhase[currChan]), BCS_LCA_d);
		repNoFlip[i] = codeChips[sampIdx];

		// For this channel, check if the nav bit flips in this sample set
		if ((idxNextNav[currChan] > 0) && (idxNextNav[currChan] < BCS_S_d)) {
			// If it does, generate the flipped replica
			if (currSamp >= idxNextNav[currChan]) {
				repFlipped[i].x = -repNoFlip[i].x;
				repFlipped[i].y = -repNoFlip[i].y;
			} else {
				repFlipped[i].x = repNoFlip[i].x;
				repFlipped[i].y = repNoFlip[i].y;
			}

		}
		// If it does not, zero the flipped array so it has no magnitude to interfere with no_flip
		else {
			repFlipped[i].x = 0;
			repFlipped[i].y = 0;
		}

		// Advance through the samples by the number of threads
		i += stride;
	}
}


/** \brief Perform element-wise multiplication of the reference array by every array in the concatinated striding array
 *
 * \param stridingArr		The single array which is applied to every consituent array in refArr
 * \param stridingArrLen	The number of elements in stridingArr
 * \param refArr			An array containing multiple smaller arrays, each of the same length, concatinated in device memory
 * \param refArrLen			The number of elements in refArr
 * \param outArr			Output: The resultant array (the same length as stridingArr)
 * \param aux				(Debug)
 *
 */
__global__ void BCS_BatchMultiply(cufftDoubleComplex *stridingArr,
		int stridingArrLen, cufftDoubleComplex *refArr, int refArrLen,
		cufftDoubleComplex *outArr, int aux) {
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int stride = blockDim.x * gridDim.x;

	int currSamp; // The index of the sample in the sequence

	// Iterate through every sample in the reference array
	while (i < refArrLen) {

		// Update index tracker
		currSamp = POSMOD(i, stridingArrLen);

		// Complex multiply these samples
		// NOTE: This gives different values than what PyGNSS does,
		// but cuCmul seems to give more accurate results than Python's *
		outArr[i] = cuCmul(stridingArr[currSamp], refArr[i]);

		// Advance through the samples by the number of threads
		i += stride;
	}
}


/** \brief Performs batch array multiplication, but selects which refArr to use based on correlation score
 *
 * \param stridingArr		An array containing multiple smaller arrays, each of the same length, concatinated in device memory
 * \param stridingArrLen	The number of elements in stridingArr
 * \param numChan			The number of arrays to process in stridingArr
 * \param refArr1			A single array used in the batch multiplication
 * \param refArr2			A single array used in the batch multiplication
 * \param refArr1Larger		Flags indicating whether to multiply by refArr1 or refArr2 based on correlation score
 * \param outArr			Output: The resultant array (the same length as stridingArr)
 * \param outArrBatchLen	Allows for padding -- the desired size of each array in the concatinated output
 *
 */
__global__ void BCS_ChoosyBatchMultiplyAndPad(cufftDoubleComplex *stridingArr,
		int stridingArrLen, int numChan, cufftDoubleComplex *refArr1,
		cufftDoubleComplex *refArr2, bool *refArr1Larger,
		int outArrBatchLen, cufftDoubleComplex *outArr) {
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int stride = blockDim.x * gridDim.x;

	int currChan; // The channel this thread is processing
	int currSamp; // The index of the sample in the sequence


	while (i < (stridingArrLen * numChan)) {

		// Update index trackers
		currChan = i / stridingArrLen;
		currSamp = POSMOD(i, stridingArrLen);

		// Choose the refArr with the higher correlation score at 0 offset for this channel
		if (refArr1Larger[currChan]) {
			// The first array scores better, mix with that
			outArr[(currChan*outArrBatchLen) + currSamp] = cuCmul(stridingArr[(currChan * stridingArrLen) + currSamp],
					refArr1[(currChan * stridingArrLen) + currSamp]);
		} else {
			// The second array scores better, mix with that
			outArr[(currChan*outArrBatchLen) + currSamp] = cuCmul(stridingArr[(currChan * stridingArrLen) + currSamp],
					refArr2[(currChan * stridingArrLen) + currSamp]);
		}

		i += stride;
	}
}






/** \brief Subtract the DC offset from the raw samples
 *
 * The DC offset must be multiplied by the Doppler wipeoff, as the input samples are intended to have had the Doppler effect removed already
 *
 * \param refArr			Pointer to an array with Doppler-wiped samples for each channel
 * \param refArrLen			The number of elements in refArr
 * \param dcOffset			The DC offset of the raw samples
 * \param dcOffsetScaleArr	The Doppler wipeoff values for each sample in each channel
 * \param destArr			Output: The DC-subtracted, Doppler-wiped samples
 *
 */
__global__ void BCS_SubtractDCOffset(cufftDoubleComplex *refArr, int refArrLen,
		cufftDoubleComplex dcOffset, cufftDoubleComplex *dcOffsetScaleArr,
		cufftDoubleComplex *destArr) {
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int stride = blockDim.x * gridDim.x;

	// Iterate through every sample in the reference array
	while (i < refArrLen) {

		// Subtract off the mean of that channel (scaling it the same way that the ref signal is scaled
		destArr[i] = cuCsub(refArr[i], cuCmul(dcOffset, dcOffsetScaleArr[i]));

		// Advance through the samples by the number of threads
		i += stride;
	}
}



/** \brief Decide between the Flipped and No-Flipped replica scores for each channel, putting the chosen scores into a composite array
 *
 * \param refArr1		The first form of scores for each channel
 * \param refArr2		The second form of scores for each channel
 * \param numChan		The number of arrays to process in stridingArr
 * \param idxNext		The index of the next navigation bit in the samples for each channel
 * \param betterScore	Output: The composite array of the chosen scores
 * \param refArr1Larger	Output: Flags indicating whether refArr1 or refArr2 was chosen
 *
 */
__global__ void BCS_ChooseCodeCorr(cufftDoubleComplex *refArr1,
		cufftDoubleComplex *refArr2, int numChan, int *idxNext,
		cufftDoubleComplex *betterScore, bool *refArr1Larger) {
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int stride = blockDim.x * gridDim.x;

	int currChan; // The channel this thread is processing


	// First, find out for every channel if the flip or noflip data should be used
	while (i < numChan) {
		// If there was no possible flip this sample set or refArr1 has the higher correlation score at 0
		// for this channel, choose refArr1 as the better score
		if (!((idxNext[i] > 0) && (idxNext[i] < BCS_S_d)) || cuCabs(refArr1[i * BCS_S_d]) > cuCabs(refArr2[i * BCS_S_d])) {
			refArr1Larger[i] = 1;
		} else { // By DeMorgan's, choose refArr2 as the better score
			refArr1Larger[i] = 0;
		}

		i += stride;
	}

	// Reset thread index so the next segment will start clean
	i = blockIdx.x * blockDim.x + threadIdx.x;

	__syncthreads();

	// Iterate through every sample in the reference array
	while (i < (numChan * BCS_S_d)) {

		// Update index tracker
		currChan = i / BCS_S_d;

		// Copy the better correlation score to the output array
		if (refArr1Larger[currChan]) {
			betterScore[i] = refArr1[i];
		} else {
			betterScore[i] = refArr2[i];
		}

		// Advance through the samples by the number of threads
		i += stride;
	}

}

/** \brief Shift the 0-frequency component of the FFT to the center of the spectrum for each array in the composite array
 * 
 * Adapted from: https://github.com/marwan-abdellah/cufftShift/blob/master/Src/CUDA/Kernels/in-place/cufftShift_1D_IP.cu
 * 
 * \param arr		The array of arrays to shift (in place, each sub-array individually)
 * \param arrLen	The number of elements in arr
 * \param batchLen	The number of elements of each sub-array in arr
 * 
 */
__global__ void BCS_cufftBatchShift(cufftDoubleComplex *arr, int arrLen,
		int batchLen) {
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int stride = blockDim.x * gridDim.x;

	int currSamp; // The index of the sample in the sequence
	int currChan; // The channel this thread is processing
	int halfBatch = batchLen / 2; // i should only track half way through each batch
	int halfArr = arrLen / 2; // Thus, i only tracks half way through the array

	// Iterate through every sample in the reference array
	while (i < halfArr) {

		// Update index tracker
		currSamp = POSMOD(i, halfBatch);
		currChan = i / halfBatch;

		// Hold the current value
		cufftDoubleComplex temp = arr[(currChan * batchLen) + currSamp];

		// Note: currChan*batchLen + currSamp == i
		// Swap the first element
		arr[(currChan * batchLen) + currSamp] = arr[(currChan * batchLen)
				+ currSamp + halfBatch];
		// Swap the second element
		arr[(currChan * batchLen) + currSamp + halfBatch] = temp;

		// Advance through the array by the number of threads
		i += stride;
	}
}


/** \brief Normalize FFT scores for each array in the composite array
 *
 * \param arr			The array of arrays to normalize (in place, each sub-array individually)
 * \param batchLen		The number of elements per subarray in arr
 * \param numBatches	The number of subarrays in arr
 *
 */
__global__ void BCS_NormalizeFFT(cufftDoubleComplex *arr, int batchLen,
		int numBatches) {

	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int stride = blockDim.x * gridDim.x;

	while (i < batchLen * numBatches) {
		arr[i].x = arr[i].x / ((float) batchLen);
		arr[i].y = arr[i].y / ((float) batchLen);

		i += stride;
	}
}





/** \brief Debug function to make viewing device memory easier. Placing a breakpoint in a kernel doesn't guarantee thread safety.
 * 		   This function makes no modifications to data, so run it on a stream to view device memory at that point.
 *
 * \param *arr		The cufftDoubleComplex array you want to view
 * \param arrLen	The size of the array you want to view
 * \param aux		Give this a unique value to identify where you are in the stream
 *
 */
__global__ void BCS_valInspector(cufftDoubleComplex *arr, int arrLen, int aux) {
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int stride = blockDim.x * gridDim.x;

	//int currChan;
	//int currSamp;

	int test;

	while (i < arrLen) {

		// You can put a breakpoint here to inspect the device arrays
		if (i == 0) {
			test = 1;
		}

		i += stride;
	}

}


/** \brief Debug function to make viewing device memory easier. Placing a breakpoint in a kernel doesn't guarantee thread safety.
 * 		   This function makes no modifications to data, so run it on a stream to view device memory at that point.
 *
 * \param *arr		The double array you want to view
 * \param arrLen	The size of the array you want to view
 * \param aux		Give this a unique value to identify where you are in the stream
 *
 */
__global__ void BCS_valInspector2(double *arr, int arrLen, int calledby) {
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int stride = blockDim.x * gridDim.x;

	//int currChan;
	//int currSamp;

	int test;

	while (i < arrLen) {

		// You can put a breakpoint here to inspect the device arrays
		if (i == 0) {
			test = 1;
		}

		i += stride;
	}

}




dsp::BatchCorrScores::BatchCorrScores() {
	ModuleName = "BatchCorrScores";
	AllocateInputs(12);
	AllocateOutputs(3);

	Started = 0;

	ConfigExpectedInput(0, "Samples", UNDEFINED_t, VALUE_CMPX, VECTORLENGTH_ANY);
	ConfigExpectedInput(1, "ValidPRNs", CHAR_t, VALUE, VECTORLENGTH_ANY);
	ConfigExpectedInput(2, "CodePhaseStart", DOUBLE_t, VALUE, VECTORLENGTH_ANY);
	ConfigExpectedInput(3, "CarrierPhaseStart", DOUBLE_t, VALUE, VECTORLENGTH_ANY);
	ConfigExpectedInput(4, "CodeFrequency", DOUBLE_t, FREQUENCY_HZ, VECTORLENGTH_ANY);
	ConfigExpectedInput(5, "CarrierFrequency", DOUBLE_t, FREQUENCY_HZ, VECTORLENGTH_ANY);
	ConfigExpectedInput(6, "cpElapsedStart", INT_t, VALUE, VECTORLENGTH_ANY);
	ConfigExpectedInput(7, "cpReference", INT_t, VALUE, VECTORLENGTH_ANY);
	ConfigExpectedInput(8, "DopplerSign", INT_t, VALUE, VECTORLENGTH_ANY);
	ConfigExpectedInput(9, "SamplingFrequency", DOUBLE_t, FREQUENCY_HZ, 1);
	ConfigExpectedInput(10, "SampleLength", DOUBLE_t, VALUE, 1);
	// Note: Phase and cpEla params are referenced to the front of the sample set!!!!

	std::clog << "[" << ModuleName << "] Configured inputs" << std::endl;

	ConfigOutput(0, "CodeScores", UNDEFINED_t, VALUE_CMPX, CUDA_DEVICE, VECTORLENGTH_ANY, NULL, 0);
	ConfigOutput(1, "CarrScores", UNDEFINED_t, VALUE_CMPX, CUDA_DEVICE,	VECTORLENGTH_ANY, NULL, 0);
	ConfigOutput(2, "NumFFTPoints", INT_t, VALUE, HOST, 1, NULL, 0);

}

dsp::BatchCorrScores::~BatchCorrScores() {
	if (Started)
		Stop();
	delete[] expectedInputs;
	delete[] inputs;
	delete[] outputs;
}

int dsp::BatchCorrScores::Start(void* cuFlowStream) {
	if (Started) {
		std::clog << "[" << ModuleName << "] Start: Already Started."
				<< std::endl;
		return 0;
	}
	std::clog << "[" << ModuleName << "] Starting ... " << std::flush;


	// Set the CUDA Stream for the GPU operations
	cuStream = (cudaStream_t*) cuFlowStream;

	// Set the initial value used in thrust reduction to 0
	signalBiasInit.x = 0;
	signalBiasInit.y = 0;

	// Set up streams for processing
	cuCheck(cudaStreamCreate(&codeStream));
	cuCheck(cudaStreamCreate(&carrStream));
	cuCheck(cudaStreamCreate(&interStream));
	cuCheck(cudaStreamCreate(&satStream));

	// Set up events for processing
	cuCheck(cudaEventCreate(&eventA));
	cuCheck(cudaEventCreate(&eventB));
	cuCheck(cudaEventCreate(&eventC));
	cuCheck(cudaEventCreate(&eventD));
	cuCheck(cudaEventCreate(&eventE));
	cuCheck(cudaEventCreate(&eventF));
	cuCheck(cudaEventCreate(&eventG));
	cuCheck(cudaEventCreate(&eventH));
	cuCheck(cudaEventCreate(&eventI));
	cuCheck(cudaEventCreate(&eventK));
	cuCheck(cudaEventCreate(&eventM));
	cuCheck(cudaEventCreate(&eventN));
	cuCheck(cudaEventCreate(&eventO));

	// Generate CA code for all PRNs
	cuCheckMSt(cudaMalloc((void**)&chipsCACode_d, sizeof(cufftDoubleComplex) * CONST_PRN_MAX * 1023));
	BCS_GenCACode<<<1, auxil::roundUpToNextPowerOfTwo(CONST_PRN_MAX), 0, *cuStream>>>(chipsCACode_d, CONST_PRN_MAX);

	// Number of samples in sample set
	S = inputs[0]->VectorLength;
	// Sampling frequency
	fs = *((double*) (inputs[9]->Data));
	// Length of one sample set
	T = *((double*) (inputs[10]->Data));
	// Round to the nearest 1us to ensure there is no trailing garbage precision
	T = round(T * 1.0e6) / 1.0e6;

	N = S * (T / 0.001); // N is number of 1ms C/A code periods in the sample set
	carrSTot = auxil::roundUpToNextPowerOfTwo(S) * 8; // Pad the carrier fft by zero-extending to 8x nearest power of 2
	truncLen = 16384;
	carrSTotTrunc = truncLen * 8; // Truncated fft: largest power of 2 within 25000 samples with 8x zero-padding for resolution




	// Copy processing parameters to device constant memory
	cuCheckMSt(cudaMemcpyToSymbol(BCS_LCA_d, &numCA, sizeof(int), 0, cudaMemcpyHostToDevice));
	cuCheckMSt(cudaMemcpyToSymbol(BCS_FS_d, &fs, sizeof(double), 0, cudaMemcpyHostToDevice));
	cuCheckMSt(cudaMemcpyToSymbol(BCS_S_d, &S, sizeof(int), 0, cudaMemcpyHostToDevice));
	cuCheckMSt(cudaMemcpyToSymbol(BCS_T_d, &T, sizeof(double), 0, cudaMemcpyHostToDevice));

	// Compute the time indices of each sample
	cuCheckMSt(cudaMalloc((void** )&timeIdcs_d, sizeof(double) * S));
	BCS_GenTimeIdcs<<<2, 512, 0, *cuStream>>>(timeIdcs_d);
	cuCheckMSt(cudaPeekAtLastError());

	// Memory allocations of buffers
	size_t size = sizeof(cufftDoubleComplex) * S;
	cuCheckMSt(cudaMalloc((void** )&rawSignal_d, size));

	// For the following buffers, note that space is allocated through PRN_MAX,
	// however, only the slots up to the number of currently tracked satellites
	// is used. If memory becomes an issue, it may be better to reduce this.
	size = size * CONST_PRN_MAX;
	codeCorrSize = size; // Remember this size while we're here since these arrays must be re-zeroed
	cuCheckMSt(cudaMalloc((void** )&dopplerWipeoff_d, size));
	cuCheckMSt(cudaMalloc((void** )&rawWiped_d, size));
	cuCheckMSt(cudaMalloc((void** )&rawZeroMeanWiped_d, size));
	cuCheckMSt(cudaMalloc((void** )&rawWipedfft_d, size));
	cuCheckMSt(cudaMalloc((void** )&replicaNoFlip_d, size));
	cuCheckMSt(cudaMalloc((void** )&replicaFlipped_d, size));
	cuCheckMSt(cudaMalloc((void** )&replicaNoFlipConjfft_d, size));
	cuCheckMSt(cudaMalloc((void** )&replicaFlippedConjfft_d, size));
	cuCheckMSt(cudaMalloc((void** )&codeCorrNoFlipfft_d, size));
	cuCheckMSt(cudaMalloc((void** )&codeCorrFlippedfft_d, size));
	cuCheckMSt(cudaMalloc((void** )&codeCorrNoFlipRaw_d, size));
	cuCheckMSt(cudaMalloc((void** )&codeCorrFlippedRaw_d, size));
	cuCheckMSt(cudaMalloc((void** )&codeCorrOut_d, size));
	cuCheckMSt(cudaMalloc((void**)&replicaNoFlipDouble_d, sizeof(double)*S*CONST_PRN_MAX));
	cuCheckMSt(cudaMalloc((void**)&replicaFlippedDouble_d, sizeof(double)*S*CONST_PRN_MAX));

	size = sizeof(int) * CONST_PRN_MAX;
	cuCheckMSt(cudaMalloc((void** )&idxNextNavBit_d, size));
	cuCheckMSt(cudaMalloc((void** )&cpCompleted_d, size));

	size = sizeof(bool) * CONST_PRN_MAX;
	cuCheckMSt(cudaMalloc((void** )&noFlipIsLarger_d, size));
	cuCheckMSt(cudaMalloc((void** )&validPRN_d, size));

	size = sizeof(cufftDoubleComplex) * CONST_PRN_MAX * carrSTot;
	cuCheckMSt(cudaMalloc((void** )&carrBaseband_d, size));
	cuCheckMSt(cudaMalloc((void** )&carrfftOut_d, size));

	size = sizeof(cufftDoubleComplex) * CONST_PRN_MAX * truncLen;
	cuCheckMSt(cudaMalloc((void**)&rawSignalTrunc_d, size));

	size = sizeof(cufftDoubleComplex) * CONST_PRN_MAX * carrSTotTrunc;
	cuCheckMSt(cudaMalloc((void**)&carrBasebandTrunc_d, size));
	cuCheckMSt(cudaMalloc((void**)&carrfftOutTrunc_d, size));

	size = sizeof(cufftDoubleComplex) * CONST_PRN_MAX;
	cuCheckMSt(cudaMalloc((void**)&rawMeanTrunc_d, size));



	// Zero out buffers for all intermediate steps
	size = sizeof(cufftDoubleComplex) * S;
	cuCheckMSt(cudaMemsetAsync((void* )rawSignal_d, 0, size, *cuStream));

	size = size * CONST_PRN_MAX;
	cuCheckMSt(cudaMemsetAsync((void* )dopplerWipeoff_d, 0, size, *cuStream));
	cuCheckMSt(cudaMemsetAsync((void* )rawWiped_d, 0, size, *cuStream));
	cuCheckMSt(cudaMemsetAsync((void* )rawZeroMeanWiped_d, 0, size, *cuStream));
	cuCheckMSt(cudaMemsetAsync((void* )rawWipedfft_d, 0, size, *cuStream));
	cuCheckMSt(cudaMemsetAsync((void* )replicaNoFlip_d, 0, size, *cuStream));
	cuCheckMSt(cudaMemsetAsync((void* )replicaFlipped_d, 0, size, *cuStream));
	cuCheckMSt(cudaMemsetAsync((void* )replicaNoFlipConjfft_d, 0, size,	*cuStream));
	cuCheckMSt(cudaMemsetAsync((void* )replicaFlippedConjfft_d, 0, size, *cuStream));
	cuCheckMSt(cudaMemsetAsync((void* )codeCorrNoFlipfft_d, 0, size, *cuStream));
	cuCheckMSt(cudaMemsetAsync((void* )codeCorrFlippedfft_d, 0, size, *cuStream));
	cuCheckMSt(cudaMemsetAsync((void* )codeCorrNoFlipRaw_d, 0, size, *cuStream));
	cuCheckMSt(cudaMemsetAsync((void* )codeCorrFlippedRaw_d, 0, size, *cuStream));

	size = sizeof(int) * CONST_PRN_MAX;
	cuCheckMSt(cudaMemsetAsync((void* )idxNextNavBit_d, 0, size, *cuStream));
	cuCheckMSt(cudaMemsetAsync((void* )cpCompleted_d, 0, size, *cuStream));

	size = sizeof(bool) * CONST_PRN_MAX;
	cuCheckMSt(cudaMemsetAsync((void* )noFlipIsLarger_d, 0, size, *cuStream));

	size = sizeof(cufftDoubleComplex) * CONST_PRN_MAX * carrSTot;
	cuCheckMSt(cudaMemsetAsync((void* )carrBaseband_d, 0, size, *cuStream));
	cuCheckMSt(cudaMemsetAsync((void* )carrfftOut_d, 0, size, *cuStream));

	size = sizeof(cufftDoubleComplex) * CONST_PRN_MAX * truncLen;
	cuCheckMSt(cudaMemsetAsync((void* )rawSignalTrunc_d, 0, size, *cuStream));

	size = sizeof(cufftDoubleComplex) * CONST_PRN_MAX * carrSTotTrunc;
	cuCheckMSt(cudaMemsetAsync((void* )carrBasebandTrunc_d, 0, size, *cuStream));
	cuCheckMSt(cudaMemsetAsync((void* )carrfftOutTrunc_d, 0, size, *cuStream));

	size = sizeof(cufftDoubleComplex) * CONST_PRN_MAX;
	cuCheckMSt(cudaMemsetAsync((void* )rawMeanTrunc_d, 0, size, *cuStream));


	// Now that the space is allocated, the output pointers can be assigned
	outputs[0].Data = codeCorrOut_d;
	outputs[0].VectorLength = (S * CONST_PRN_MAX) / (int) (N + 0.5);
	outputs[1].Data = carrfftOut_d;
	outputs[1].VectorLength = CONST_PRN_MAX * carrSTot;
	outputs[2].Data = (void*) &carrSTot;
	outputs[2].VectorLength = 1;

	// Make sure all GPU tasks have completed before continuing
	cuCheckMSt(cudaDeviceSynchronize());
	cuCheckMSt(cudaStreamSynchronize(*cuStream));

	// Signifies that the next call to update() will be the first after start()
	Started = 1;
	FirstUpdate = 1;

	std::clog << "Started." << std::endl;
	return 0;
}

int dsp::BatchCorrScores::Stop() {
	int ret = 0;

	if (Started == 0) {
		std::clog << "[" << ModuleName << "] Stop: Wasn't Started."
				<< std::endl;
		return 0;
	}
	std::clog << "[" << ModuleName << "] Stopping ... " << std::flush;

	// Destroy processing streams
	cuCheckV(cudaStreamDestroy(codeStream));
	cuCheckV(cudaStreamDestroy(carrStream));
	cuCheckV(cudaStreamDestroy(interStream));
	cuCheckV(cudaStreamDestroy(satStream));

	// Destroy processing events
	cuCheckV(cudaEventDestroy(eventA));
	cuCheckV(cudaEventDestroy(eventB));
	cuCheckV(cudaEventDestroy(eventC));
	cuCheckV(cudaEventDestroy(eventD));
	cuCheckV(cudaEventDestroy(eventE));
	cuCheckV(cudaEventDestroy(eventF));
	cuCheckV(cudaEventDestroy(eventG));
	cuCheckV(cudaEventDestroy(eventH));
	cuCheckV(cudaEventDestroy(eventI));
	cuCheckV(cudaEventDestroy(eventK));
	cuCheckV(cudaEventDestroy(eventM));
	cuCheckV(cudaEventDestroy(eventN));
	cuCheckV(cudaEventDestroy(eventO));

	// Free device memory
	cuCheckMSp(cudaFree((void* )rawSignal_d));
	cuCheckMSp(cudaFree((void* )PRNs_d));
	cuCheckMSp(cudaFree((void* )codePhase_d));
	cuCheckMSp(cudaFree((void* )carrierPhase_d));
	cuCheckMSp(cudaFree((void* )codeFrequency_d));
	cuCheckMSp(cudaFree((void* )carrierFrequency_d));
	cuCheckMSp(cudaFree((void* )cpElapsed_d));
	cuCheckMSp(cudaFree((void* )cpReference_d));
	cuCheckMSp(cudaFree((void* )timeIdcs_d));
	cuCheckMSp(cudaFree((void* )dopplerWipeoff_d));
	cuCheckMSp(cudaFree((void* )rawWiped_d));
	cuCheckMSp(cudaFree((void* )rawZeroMeanWiped_d));
	cuCheckMSp(cudaFree((void* )rawWipedfft_d));
	cuCheckMSp(cudaFree((void* )carrBaseband_d));
	cuCheckMSp(cudaFree((void* )replicaNoFlip_d));
	cuCheckMSp(cudaFree((void* )replicaFlipped_d));
	cuCheckMSp(cudaFree((void* )replicaNoFlipConjfft_d));
	cuCheckMSp(cudaFree((void* )replicaFlippedConjfft_d));
	cuCheckMSp(cudaFree((void* )idxNextNavBit_d));
	cuCheckMSp(cudaFree((void* )cpCompleted_d));
	cuCheckMSp(cudaFree((void* )replicaNoFlipConjfft_d));
	cuCheckMSp(cudaFree((void* )replicaFlippedConjfft_d));
	cuCheckMSp(cudaFree((void* )codeCorrNoFlipfft_d));
	cuCheckMSp(cudaFree((void* )codeCorrFlippedfft_d));
	cuCheckMSp(cudaFree((void* )codeCorrNoFlipRaw_d));
	cuCheckMSp(cudaFree((void* )codeCorrFlippedRaw_d));
	cuCheckMSp(cudaFree((void* )codeCorrNoFlipPerCP_d));
	cuCheckMSp(cudaFree((void* )codeCorrFlippedPerCP_d));
	cuCheckMSp(cudaFree((void* )noFlipIsLarger_d));
	cuCheckMSp(cudaFree((void* )codeCorrOut_d));
	cuCheckMSp(cudaFree((void* )carrfftOut_d));
	cuCheckMSp(cudaFree((void* )validPRN_d));
	cuCheckMSp(cudaFree((void* )txTime_d));
	cuCheckMSp(cudaFree((void* )TOWcpReference_d));
	cuCheckMSp(cudaFree((void* )replicaNoFlipDouble_d));
	cuCheckMSp(cudaFree((void* )replicaFlippedDouble_d));
	cuCheckMSp(cudaFree((void* )currSatStates_d));
	cuCheckMSp(cudaFree((void* )rawSignalTrunc_d));
	cuCheckMSp(cudaFree((void* )carrBasebandTrunc_d));
	cuCheckMSp(cudaFree((void* )carrfftOutTrunc_d));
	cuCheckMSp(cudaFree((void* )rawMeanTrunc_d));



	// Remove fft plans
	cufftCheckMSp(cufftDestroy(codeStreamfftHandle));
	cufftCheckMSp(cufftDestroy(carrStreamfftHandle));
	cufftCheckMSp(cufftDestroy(interStreamfftHandle));

	Started = 0;
	std::clog << "Stopped." << std::endl;

	return ret;
}

int dsp::BatchCorrScores::Update(void* cuFlowStream) {
	if (Started == 0) {
		std::cerr << "[" << ModuleName
				<< "] Error: Update() Failed due to batch correlator not initialized"
				<< std::endl;
		return -1;
	}

	if (syncFailed) {
		std::cerr << "[" << ModuleName
				<< "] Error: Update() Failed due to stream synchronization error"
				<< std::endl;
		return -1;
	}

	// Set the CUDA Stream for the fft GPU operations
	if (FirstUpdate) {

		// Determine how many channels are being tracked (cuChanMgr outputs not set in BCS.Start())
		numChan = inputs[1]->VectorLength; // number of PRNs currently being tracked

		// Get pointers to the channel parameters (cuChanMgr outputs not set in BCS.Start())
		codePhase_d 		= (double*)(inputs[2]->Data);
		carrierPhase_d 		= (double*)(inputs[3]->Data);
		codeFrequency_d 	= (double*)(inputs[4]->Data);
		carrierFrequency_d 	= (double*)(inputs[5]->Data);
		cpElapsed_d 		= (int*)(inputs[6]->Data);
		cpReference_d 		= (int*)(inputs[7]->Data);
		dopplerSign_d 		= (int*)(inputs[8]->Data);
		PRNs_d 				= (uint8_t*)(inputs[1]->Data);


		// Create plan for the replica correlation
		int rank = 1;           // Dimensionality of transform (1, 2, or 3)
		int n[] = { S };          // Array of size rank, the size of each dimension
		int istride = 1;        // Distance between successive input elements
		int ostride = 1;        // Distance between successive output elements
		int idist = S; // Distance between first element of two consecutive batches in input data
		int odist = S; // Distance between first element of two consecutive batches in output data
		int *inembed = NULL; // Pointer of size rank that indicates storage dimension of input data
		int *onembed = NULL; // Pointer of size rank that indicates storage dimension of output data
		int batch = numChan; // Batch size for this transform

		// Make one handle for replica fft and one for raw fft (since they'll be on different threads)
		cufftCheckMSt(cufftPlanMany(&codeStreamfftHandle, rank, n, inembed, istride,
						idist, onembed, ostride, odist, CUFFT_Z2Z, batch));
		cufftCheckMSt(cufftPlanMany(&interStreamfftHandle, rank, n, inembed, istride,
						idist, onembed, ostride, odist, CUFFT_Z2Z, batch));
		// Resize for the carrfft
		n[0] = {carrSTot};
		idist = carrSTot;
		odist = carrSTot;
		cufftCheckMSt(cufftPlanMany(&carrStreamfftHandle, rank, n, inembed, istride,
						idist, onembed, ostride, odist, CUFFT_Z2Z, batch));
		// NOTE: Consider cufftMakePlanMany if fs and/or S/N/T need to be changed during Update()!
		//       Current architecture keeps these values fixed by executive design decision.

		// Set the fft streams
		cufftCheckM(cufftSetStream(codeStreamfftHandle, codeStream));
		cufftCheckM(cufftSetStream(interStreamfftHandle, interStream));
		cufftCheckM(cufftSetStream(carrStreamfftHandle, carrStream));


		FirstUpdate = 0;
	}


	// Get the input data accessible
	rawFixed = (int16_t*) inputs[0]->Data; // Buffer pointer needs updated every iteration, since samples are being read in to a set of buffers in the background


	// eventN: Compute the location of the next nav bit
	// <<<   >>> -- needs one thread per PRN -- one block of 64 threads for 37 chans
	BCS_NavBitBoundary<<<1, auxil::roundUpToNextPowerOfTwo(CONST_PRN_MAX), 0,
			codeStream>>>(cpElapsed_d, cpReference_d, codePhase_d,
			codeFrequency_d, numChan, idxNextNavBit_d);
	//BCS_valInspector2<<<2, 1024, 0, codeStream>>>(validPRN_d, 37, 1); // Debug
	cudaEventRecord(eventN, codeStream);


	// eventA: Convert the samples into cufftDoubleComplex type
	// Threads of 64 to match NavBitBoundary, (63+1)x64 = 4096 to fill TX2 thread count
	BCS_Load<<<64, 64, 0, *cuStream>>>(rawFixed, rawSignal_d, S);
	cudaEventRecord(eventA, *cuStream);
	//BCS_valInspector<<<1, 128, 0, *cuStream>>>(rawSignal_d, 37, 100); // Debug


	// eventB: Calculate the average value of the sample set to remove the DC offset
	// (DC bias gives a large response at 0Hz in fft)
	cudaStreamWaitEvent(carrStream, eventA, 0);
	rawMean = ComplexDivide(thrust::reduce(thrust::cuda::par.on(carrStream), rawSignal_d,
					rawSignal_d + S, signalBiasInit, addCufftComplex()), (float) S);
	cudaEventRecord(eventB, carrStream);


	// eventC: Compute the Doppler wipeoff value given the channel parameters
	BCS_ComputeDopplerWipeoff<<<8, 256, 0, interStream>>>(carrierFrequency_d,
			carrierPhase_d, S, numChan, timeIdcs_d, dopplerWipeoff_d);
	//BCS_valInspector<<<2, 1024, 0, interStream>>>(dopplerWipeoff_d, 37, 110); // Debug
	cudaEventRecord(eventC, interStream);


	// eventD: Create the replica signal given the channel parameters
	// Note: The replica signal is created using PREVIOUS ITERATION's phase channel parameter values!
	//       This is because time is always referenced to the END of the current sample set --
	//       because we have the full current sample set at this moment, so we're at the end of that set.
	//       However, the replica is generated by counting FORWARD by the way the algorithm is written.
	//       Since the channel parameters TimeUpdate doesn't update frequencies, the MeasUpdate channel params
	//       may be used here, counting forward from the code phase at the END of the PREVIOUS iteration.
	//       So, even though the kernel arguments are (k-1|k-1), counting forward from (k-1|k-1) is equivalent
	//       to counting backwards from (k|k-1) when generating the replicas!
	//       However, the channel parameters should be advanced to (k|k-1) immediately after this (eventO).
	cudaStreamWaitEvent(codeStream, eventN, 0);
	BCS_ComputeCodeReplica<<<8, 256, 0, codeStream>>>(codeFrequency_d,
				codePhase_d, S, numChan, PRNs_d, timeIdcs_d, idxNextNavBit_d,
				chipsCACode_d, replicaNoFlip_d, replicaFlipped_d);
	//BCS_valInspector<<<2, 1024, 0, codeStream>>>(replicaNoFlip_d, 37, 121); // Debug
	//BCS_valInspector<<<2, 1024, 0, codeStream>>>(replicaFlipped_d, 37, 122); // Debug
	cudaEventRecord(eventD, codeStream);


	// eventE: Compute the fftconj of the replica signal (may be expensive since 2x prn*VectorLength ffts)
	cudaStreamWaitEvent(codeStream, eventD, 0);
	// conjfft of no_flip replica
	cufftCheckMSt(cufftExecZ2Z(codeStreamfftHandle, replicaNoFlip_d, replicaNoFlipConjfft_d, CUFFT_FORWARD));
	cufftCheckMSt(cufftExecZ2Z(codeStreamfftHandle, replicaFlipped_d, replicaFlippedConjfft_d, CUFFT_FORWARD));
	//BCS_valInspector<<<2, 1024, 0, codeStream>>>(replicaNoFlipConjfft_d, 37, 131); // Debug
	BCS_ComplexVectorTranspose<<<8, 256, 0, codeStream>>>(
			replicaNoFlipConjfft_d, S * numChan);
	BCS_ComplexVectorTranspose<<<8, 256, 0, codeStream>>>(
			replicaFlippedConjfft_d, S * numChan);
	//BCS_valInspector<<<2, 1024, 0, codeStream>>>(replicaFlippedConjfft_d, 37, 134); // Debug
	cudaEventRecord(eventE, codeStream);


	// eventK: Apply Doppler wipeoff to raw
	cudaStreamWaitEvent(interStream, eventA, 0);
	cudaStreamWaitEvent(interStream, eventC, 0);
	BCS_BatchMultiply<<<8, 256, 0, interStream>>>(rawSignal_d, S,
			dopplerWipeoff_d, S * numChan, rawWiped_d, 0);
	//BCS_valInspector<<<2, 1024, 0, interStream>>>(rawWiped_d, 37, 141); // Debug
	cudaEventRecord(eventK, interStream);

	// eventF: Compute fft of raw signal
	cudaStreamWaitEvent(interStream, eventK, 0);
	cufftCheckMSt(cufftExecZ2Z(interStreamfftHandle, rawWiped_d, rawWipedfft_d, CUFFT_FORWARD));
	//BCS_valInspector<<<2, 1024, 0, interStream>>>(rawWipedfft_d, 64, 7); // Debug
	cudaEventRecord(eventF, interStream);


	// eventG: Compute fft(raw)*fftconj(replica)
	cudaStreamWaitEvent(codeStream, eventE, 0);
	cudaStreamWaitEvent(codeStream, eventF, 0);
	BCS_BatchMultiply<<<8, 256, 0, codeStream>>>(rawWipedfft_d, S * numChan,
			replicaNoFlipConjfft_d, S * numChan, codeCorrNoFlipfft_d, 1);
	BCS_BatchMultiply<<<8, 256, 0, codeStream>>>(rawWipedfft_d, S * numChan,
			replicaFlippedConjfft_d, S * numChan, codeCorrFlippedfft_d, 2);
	//BCS_valInspector<<<2, 1024, 0, codeStream>>>(codeCorrNoFlipfft_d, 37, 161); // Debug
	//BCS_valInspector<<<2, 1024, 0, codeStream>>>(codeCorrFlippedfft_d, 37, 162); // Debug
	cudaEventRecord(eventG, codeStream);


	// ***eventH***: Compute code-corr by ifft
	cudaStreamWaitEvent(codeStream, eventG, 0);
	// Compute the ifft
	cufftCheckMSt(cufftExecZ2Z(codeStreamfftHandle, codeCorrNoFlipfft_d, codeCorrNoFlipRaw_d, CUFFT_INVERSE));
	cufftCheckMSt(cufftExecZ2Z(codeStreamfftHandle, codeCorrFlippedfft_d, codeCorrFlippedRaw_d, CUFFT_INVERSE));
	// Normalize IFFT to match PyGNSS (shouldn't be needed outside of debug/verification, but done for completeness sake)
	BCS_NormalizeFFT<<<8, 256, 0, codeStream>>>(codeCorrNoFlipRaw_d, S, numChan);
	BCS_NormalizeFFT<<<8, 256, 0, codeStream>>>(codeCorrFlippedRaw_d, S, numChan);
	// NOTE: the following is different from PyGNSS!
	// Instead of averaging all values with the same rc value, this algorithm will keep all sample points.
	// This will allow an error greater than 300m in the direction of a given satellite for the candidate grid positions.
	// To make the following function like PyGNSS, BCS_AvgCorrespondingSamples must be reworked so that scores corresponding
	// to the same code phase shift are averaged. Then, the remaining kernels must accomodate the correct size. Also, BCM would need updated.
	BCS_ChooseCodeCorr<<<8, 256, 0, codeStream>>>(codeCorrNoFlipRaw_d,
			codeCorrFlippedRaw_d, numChan, idxNextNavBit_d, codeCorrOut_d,
			noFlipIsLarger_d);
	BCS_cufftBatchShift<<<8, 256, 0, codeStream>>>(codeCorrOut_d, numChan * S, S);
	//BCS_valInspector<<<2, 1024, 0, codeStream>>>(codeCorrOut_d, 64, 16); // Debug
	cudaEventRecord(eventH, codeStream);


	// eventM: Compute (raw*DopplerWipeoff) - (mean*DopplerWipeoff)
	// Old way
	cudaStreamWaitEvent(carrStream, eventK, 0);
	BCS_SubtractDCOffset<<<8, 256, 0, carrStream>>>(rawWiped_d, numChan * S,
			rawMean, dopplerWipeoff_d, rawZeroMeanWiped_d);
	cudaEventRecord(eventM, carrStream);


	// eventI: Compute (raw - mean(raw)) * replica * DopplerWipeoff
	cudaStreamWaitEvent(carrStream, eventD, 0);
	cudaStreamWaitEvent(carrStream, eventM, 0);
	cudaStreamWaitEvent(carrStream, eventH, 0); // Wait for codeCorr to determine if flip is needed
	BCS_ChoosyBatchMultiplyAndPad<<<8, 256, 0, carrStream>>>(rawZeroMeanWiped_d,
			S, numChan, replicaNoFlip_d, replicaFlipped_d, noFlipIsLarger_d,
			carrSTot, carrBaseband_d);
	//BCS_valInspector<<<2, 1024, 0, carrStream>>>(carrBaseband_d, 37, 191); // Debug
	cudaEventRecord(eventI, carrStream);


	// ***eventJ***: Compute carr-fft by fft
	cudaStreamWaitEvent(carrStream, eventI, 0);
	cufftCheckMSt(cufftExecZ2Z(carrStreamfftHandle, carrBaseband_d, carrfftOut_d, CUFFT_FORWARD));
	BCS_cufftBatchShift<<<8, 256, 0, carrStream>>>(carrfftOut_d, numChan * carrSTot, carrSTot);
	//BCS_valInspector<<<2, 1024, 0, carrStream>>>(carrfftOut_d, 37, 202); // Debug




	// Wait for all streams to complete their tasks (so no processing gets overrun next iteration)
	// Note: This blocks the host thread until completion, and, since modules
	// are executed sequentially in a flow, these buffers are guaranteed
	// to have good new data in them once clearing these gates
	// Note note: DPMeas will also need to block at the end of the update
	// before continuing so this module doesn't overwrite the new values
	cuCheckMSt(cudaStreamSynchronize(codeStream));
	cuCheckMSt(cudaStreamSynchronize(carrStream));
	cuCheckMSt(cudaStreamSynchronize(interStream));
	cuCheckMSt(cudaStreamSynchronize(*cuStream));




	if (syncFailed) {
		std::cerr << "[" << ModuleName
				<< "] Error: Update() Failed due to stream synchronization error"
				<< std::endl;
		return -1;
	}

	return 0;
}

cufftDoubleComplex dsp::BatchCorrScores::ComplexDivide(cufftDoubleComplex sum,
		float S) {
	cufftDoubleComplex temp;
	temp.x = sum.x / S;
	temp.y = sum.y / S;
	return temp;
}


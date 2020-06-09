
#ifndef INC__BATCHCORRMANIFOLD_H_
#define INC__BATCHCORRMANIFOLD_H_


#include <stdint.h>
#include "dsp.h"
#include "module.h"
#include "auxil.h"
#include "ephhelper.h"
#include "statehelper.h"
#include "gridhelper.h"
#include "consthelper.h"

#define BCM_NUM_GRID_STATES 4
// Update MAX_GRID_SIZE to consider 2*4D grids
#define BCM_MAX_GRID_SIZE (75*75*75*75 + 75*75*75*75)

#define BCM_SQR(x)	((x)*(x))
#define BCM_CUBE(x)	((x)*(x)*(x))

#define BCM_RTOL_KEPLER 1e-12		/* relative tolerance for Kepler equation */
#define BCM_MAX_ITER_KEPLER 10			/* max number of iterations for Kepler equations */


class dsp::BatchCorrManifold : public dsp::Module {

	public:

		BatchCorrManifold();

		~BatchCorrManifold();

		/** \brief Any other initialization before running. */
		int Start(void* cuFlowStream);

		/** \brief Any de-initialization before de-constructor. */
		int Stop(void);

		/** \brief Update function to run every flow loop */
		int Update(void* cuFlowStream);

	protected:

		/**
		 *  Functions
		 */


		/**
		 *  Class helper variables
		 */
		char Started = 0;
		char FirstUpdate = 1;
		cudaStream_t *cuStream;
		cudaStream_t posStream;
		cudaStream_t velStream;



		/**
		 *  Inputs
		 */
		int stateDim = -1;
		double *txTimePtr_d;
		bool *PRNs_d;
		double *codeFrequency_d;
		double *carrFrequency_d;
		double *codePhase_d;
		double *carrPhase_d;
		int *cpRefTOW_d;
		int *cpElapsedEnd_d;
		int *cpRef_d;
		cufftDoubleComplex *codeScores_d;
		cufftDoubleComplex *carrScores_d;
		double *xCurr_d;
		int *dopplerSign_d;
		double *rxTimePtr;
		int *numfftPointsPtr;
		dsp::utils::state_t<double> *satStates_d;
		dsp::utils::state_t<double> *satStatesOld_d;
		double *enu2ecefMat_d;



		/**
		 *  Params
		 */
		int 							posGridDimSizeParam = 0;
		int 							velGridDimSizeParam = 0;
		float 							gridDimSpacingParam = 0;
		dsp::utils::ManifoldGridTypes 	gridTypeParam;
		int LPower;
		bool							loadPosGrid = 0;



		/**
		 *  Outputs
		 */
		double *zVal_d;
		double *RVal_d;


		/**
		 *  Intermediate variables
		 */
		int numChan;
		// Descriptors of grid size
		int gridSize[2];
		int gridDims[2*BCM_NUM_GRID_STATES];
		int gridHalfIdx[2*BCM_NUM_GRID_STATES];
		// Device pointer for the score of each point on the grid
		double *gridPosScores_d;
		double *gridVelScores_d;
		// Device pointer for the dX of the centerpoint on the grid
		dsp::utils::statePosManifold_t<double> *gridPosLocs_d;
		dsp::utils::stateVelManifold_t<double> *gridVelLocs_d;
		// Time grid for batch satellite calculation
		double *timeGrid_d;
		// Spacing between points on the grid
		double gridDimSpacingArr[2*BCM_NUM_GRID_STATES];
		double *gridDimSpacing_d;
		// Timing
		double *sampleLengthPtr;
		double *samplingFrequencyPtr;
		int numSamps;
		double currSamplingFreq;

		// The intermediate results from the mainfolds to be averaged by reduction sum
		dsp::utils::weightState_t<double> *weightedPosStates_d;
		dsp::utils::weightState_t<double> *weightedVelStates_d;

		// Pointers to the ML state's score
		double *gridPosMaxPtr_d;
		double *gridVelMaxPtr_d;

		//int threadsPerPosManiBlock = 64;
		//int threadsPerVelManiBlock = 64;
		//int threadsPerMeasRedBlock = 8;



		double *gridPosScores_h;

        static const unsigned int FilenameCapacity = 1024;
        char Filename[FilenameCapacity] = "";
        char loadPosGridFilename[FilenameCapacity] = "";


	    // NOTE: If you change threadsPer___ManiBlock, you need to change the shared memory allocation in the kernel itself to match!
	    //       Currently hacking a compile time-known value in the kernel so weird extern allocation stuff can be ignored.
		int threadsPerPosManiBlock = 64;
		int threadsPerVelManiBlock = 64;
	    // NOTE: If you change threadsPerMeasReadBlock, you need to change the shared memory allocation in the kernel itself to match!
	    //       Currently hacking a compile time-known value in the kernel so weird extern allocation stuff can be ignored.
		int threadsPerMeasRedBlock = 8;


};

#endif

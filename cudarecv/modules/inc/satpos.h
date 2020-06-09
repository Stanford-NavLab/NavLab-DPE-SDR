
#ifndef INC__SATPOS_H_
#define INC__SATPOS_H_


#include <stdint.h>
#include "dsp.h"
#include "module.h"
#include "auxil.h"
#include "ephhelper.h"
#include "statehelper.h"

// BatchCorrManifold expects these values to stay constant throughout each run
// (could make these dynamically sized, but was deemed not worth the time for initial release)
#define SATPOS_NUM_BATCHES 3
#define SATPOS_BATCH_MIN 0.1
#define SATPOS_BATCH_FREQ 5000
#define SATPOS_BATCH_NUM_POS ((int)(SATPOS_BATCH_MIN*60*SATPOS_BATCH_FREQ))
#define SATPOS_SHARED_SIZE (sizeof(int)*(CONST_PRN_MAX+1))

#define SATPOS_SQR(x)	((x)*(x))
#define SATPOS_CUBE(x)	((x)*(x)*(x))

#define SATPOS_RTOL_KEPLER 1e-12		/* relative tolerance for Kepler equation */
#define SATPOS_MAX_ITER_KEPLER 10			/* max number of iterations for Kepler equations */

class dsp::SatPos : public dsp::Module {

	public:

	// Holds a block of pre-calculated satellite positions and related info
	typedef struct {
		double start_rxTime_a;					// The rxTime_a of satPos_d[0]
		double pos_period_sec = 0.001;			// Time between position estimates (sec)
		dsp::utils::state_t<double>* satPos_d;	// Pointer to the positions stored on device memory

	} SatPosBlock_t;
		
		SatPos();
		
		~SatPos();

		/** \brief Any other initialization before running. */
		int Start(void* cuFlowStream);

		/** \brief Any de-initialization before de-constructor. */
		int Stop(void);

		/** \brief Update function to run every flow loop */
		int Update(void* cuFlowStream);

	protected:

		cudaStream_t *cuStream;

		/**
		 *  Functions
		 */
		int FindClosestEph(dsp::utils::navLite_t *nav, double ephTOW);


		/**
		 *  Class helper variables
		 */
		char Started = 0;


		/**
		 *  Inputs
		 */
		double initRXTime_a = -1;
		dsp::utils::ephSet_t *initNavPtr;
		//int initNavSize;


		/**
		 *  Params
		 */
		double oldEphTermTime;


		/**
		 *  Outputs
		 */
		// The batch-calculated satellite positions
		dsp::utils::satState_t *satPosAll_d;
		dsp::utils::satState_t *satPos_d[SATPOS_NUM_BATCHES]; // Don't do arrays of pointers on device; cuda doesn't like this
		dsp::utils::satState_t **satPosPtr_d;
		dsp::utils::state_t<double> *satPosPtrArrAll_d;
		//dsp::utils::state_t<double> *satPosPtrArr_d[SATPOS_NUM_BATCHES];
		//dsp::utils::state_t<double>* satLoc_d[5*60*100][3];
		int satPosArrSize;
		double satPosBatchTime;
	

		/**
		 *  Intermediate variables
		 */
		cudaStream_t satPosStream;
		int prevEphMinute;
		std::vector<dsp::utils::ephSet_t> ephSetVec;
		int ephSetSize;
		dsp::utils::ephSet_t *ephSetArr;
		dsp::utils::ephSet_t *ephSetPtr_d;
		dsp::utils::ephSet_t *debugEph_d;


		dsp::utils::navLite_t *navPtr_d;

		double currRxTime;
		double currFurthestBatchTime;
		int currOldestBatchIdx;



};

#endif

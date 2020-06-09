
#ifndef INC__CUCHANMGR_H_
#define INC__CUCHANMGR_H_


#include <stdint.h>
#include "dsp.h"
#include "module.h"
#include "auxil.h"
#include "ephhelper.h"
#include "statehelper.h"

#define CHM_RTOL_KEPLER 1e-12		/* relative tolerance for Kepler equation */
#define CHM_MAX_ITER_KEPLER 10			/* max number of iterations for Kepler equations */

#define CHM_SQR(x)	((x)*(x))
#define CHM_CUBE(x)	((x)*(x)*(x))

#define CHM_RTOL_KEPLER 1e-12		/* relative tolerance for Kepler equation */
#define CHM_MAX_ITER_KEPLER 10			/* max number of iterations for Kepler equations */

class dsp::cuChanMgr : public dsp::Module {

	public:

		cuChanMgr();

		~cuChanMgr();

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
		double rxTime = -1;
		double *xk1k1_d;
		double *xkk1_d;
		double *timeGrid_d;
		int timeGridDim;
		int numChan;
		// The length of one sample set
		double T;
		//int initNavSize;


		/**
		 *  Params
		 */
		int dopplerSign;


		/**
		 *  Outputs
		 */



		/**
		 *  Intermediate variables
		 */
        // The number of chips per cp
        const int numCA = CONST_L_CA;

		std::vector<dsp::utils::ephSet_t> ephSetVec;
		int ephSetSize;
		dsp::utils::ephSet_t *ephSetArr;
		dsp::utils::ephSet_t *ephSetPtr_d;
		double *codePhaseStart_d; 		// Referenced wrt the start of the sample set
		double *carrierPhaseStart_d;	// Referenced wrt the start of the sample set
		double *codePhaseEnd_d;			// Referenced wrt the end of the sample set
		double *carrierPhaseEnd_d;		// Referenced wrt the end of the sample set
		double *codeFrequency_d;
		double *carrierFrequency_d;
		int *cpElapsedStart_d;			// Referenced wrt the start of the sample set
		int *cpElapsedEnd_d;			// Referenced wrt the end of the sample set
		int *cpReference_d;
		int *TOWcpReference_d;
		int *dopplerSign_d;
		dsp::utils::state_t<double> *satStates_d;
		uint8_t *PRNs_d;
		double *txTime_d;

		int *ephToUse_d;

		dsp::utils::state_t<double> *batchSatStates_d;
		double *enu2ecefMat_d;

		int prepBlockCount;



	};

	#endif

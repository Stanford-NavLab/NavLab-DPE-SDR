
#ifndef INC__BATCHCORRSCORES_H_
#define INC__BATCHCORRSCORES_H_

#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <cuComplex.h>
#include <cufft.h>
#include <stdint.h>
#include "dsp.h"
#include "module.h"
//#include "constant.h"
#include "ephhelper.h"
#include "statehelper.h"

#define BCS_SQR(x)	((x)*(x))
#define BCS_CUBE(x)	((x)*(x)*(x))

#define BCS_RTOL_KEPLER 1e-12		/* relative tolerance for Kepler equation */
#define BCS_MAX_ITER_KEPLER 10			/* max number of iterations for Kepler equations */


/** \brief Batch Score Correlator module. */
class dsp::BatchCorrScores : public dsp::Module {

    public:

		BatchCorrScores();

        ~BatchCorrScores();

        /** \brief Any other initialization before running. */
        int Start(void* cuFlowStream);

        /** \brief Any de-initialization before de-constructor. */
        int Stop(void);
    
        /** \brief Update function to run every flow loop */
        int Update(void* cuFlowStream);

    protected:
    
        // Functions

        /** \brief Divide both values of a complex number */
        cufftDoubleComplex ComplexDivide(cufftDoubleComplex sum, float S);

        /** Buffer to store all 37 CA codes. */
        cufftDoubleComplex hostCACode[37*1023];




        // CUDA functionality

        /** \brief Input stream for updates */
        cudaStream_t *cuStream;

        /** \brief Streams for code corr, carr fft, and intermediate operations. */
        cudaStream_t codeStream, carrStream, interStream, satStream;
        
        /** \brief Events for the different steps (see BatchCorr parallelism for meaning). */
        cudaEvent_t eventA, eventB, eventC, eventD, eventE, eventF, eventG, eventH,
                    eventI, eventK, eventM, eventN, eventO;

        // fft plan handles
        cufftHandle codeStreamfftHandle;
        cufftHandle interStreamfftHandle;
        cufftHandle carrStreamfftHandle;
        cufftHandle codeStreamfftHandle2;



        // Local helper values

        // The number of chips per cp
        const int numCA = CONST_L_CA;
        
        // For use in thrust reduction, the initial value of the sum
        cufftDoubleComplex signalBiasInit;
                 
        /** Sampling frequency */
        double fs;
        
        /** Length in time of one sample set */
        double T;
        
        /** Number of samples in one sample set */
        int S;
        
        /** Number of code periods in this sample set */
        double N;
        
        /** Size of the carrier fft in number of samples */
        int carrSTot;
        int truncLen = 16384;
        int carrSTotTrunc;
        
        /** Size of codecorr array */
        int codeCorrSize;
        int codeCorrCount;

        /** Length of CA code. */
        int LC = 1023;




        // Parameters
		int dopplerSign;




        // Inputs access
        int16_t *rawFixed;
        int numChan;
        uint8_t *PRNs_d;
        double *codePhase_d;
        double *carrierPhase_d;
        double *codeFrequency_d;
        double *carrierFrequency_d;
        int *cpElapsed_d;
        int *cpReference_d;
        int *TOWcpReference_d;
        int *dopplerSign_d;
		dsp::utils::satState_t **satPosPtrs_d;
		double *xCurr_d;

		// Debug
		dsp::utils::ephSet_t *debugEph_d;




        // Intermediate variables

        /** Flag to signify if the call to Update() is the first call after Start()
          * Purpose is to configure a cuFFT plan to the flow cuda stream. */
        char FirstUpdate = 1;

        /** Flag to signify if module is started. */
        char Started = 0;
        
        /** Flag to signify if thread synchronization in update failed */
        char syncFailed = 0;

        /** GPS time reference for delta_t estimates */
        double rxTime;

        /** Flag to signify if the fft needs to be resized (eg, new sats acquired) */
        bool updatefft = false;

        // Format: MAX_PRN * 1023
        cufftDoubleComplex *chipsCACode_d;

        // Format: Prns * VectorLength
        //cufftDoubleComplex *digiCACode_d;

        // Format: Prns * VectorLength
        cufftDoubleComplex *digiCACodeFFTconj_d;

        // The raw samples, parsed as cufftDoubleComplex values
        // Format: VectorLength
        cufftDoubleComplex *rawSignal_d;
        
        // The raw samples, truncated for a smaller carr_fft that doesn't need the navbit
        cufftDoubleComplex *rawSignalTrunc_d;

        // Debug: host copy of the raw samples
        cufftDoubleComplex *rawSignal;

        // The time wrt the first sample of each sample
        // Format: VectorLength
        double *timeIdcs_d;
        
        // The Doppler wipeoff signal associated with each sample given channel params
        // Format: Prns * VectorLength
        cufftDoubleComplex *dopplerWipeoff_d;
        
        // Raw mixed with Doppler wipeoff
        // Format: Prns * VectorLength
        cufftDoubleComplex *rawWiped_d;
        
        // (raw - mean(raw)) mixed with Doppler wipeoff
        // Format: Prns * VectorLength
        cufftDoubleComplex * rawZeroMeanWiped_d;
        
        // mean(raw)
        cufftDoubleComplex rawMean;
        cufftDoubleComplex rawMeanTrunc[CONST_PRN_MAX];
        cufftDoubleComplex* rawMeanTrunc_d;
        size_t rawMeanTruncSize = sizeof(cufftDoubleComplex) * CONST_PRN_MAX;

        // fft of Doppler wipeoff times raw
        // Format: Prns * VectorLength
        cufftDoubleComplex *rawWipedfft_d;
        
        // (raw - mean(raw)) * DopplerWipeoff * replica
        // Format: Prns * VectorLength
        cufftDoubleComplex *carrBaseband_d;
        cufftDoubleComplex *carrBasebandTrunc_d;
        
        // The code replica that should be received given channel params
        // Format: Prns * VectorLength
        cufftDoubleComplex *replicaNoFlip_d;
        cufftDoubleComplex *replicaFlipped_d;
        double *replicaNoFlipDouble_d;
        double *replicaFlippedDouble_d;
        
        // The complex conjugate of the fft of the replica signal
        // Format: Prns * VectorLength
        cufftDoubleComplex *replicaNoFlipConjfft_d;
        cufftDoubleComplex *replicaFlippedConjfft_d;
        
        // The fft of the code-correlated signal for every evaluated code shift
        // Format: Prns * VectorLength
        cufftDoubleComplex *codeCorrNoFlipfft_d;
        cufftDoubleComplex *codeCorrFlippedfft_d;
        
        // The scores of the code-correlated signal for every evaluated code shift
        // Format: Prns * VectorLength
        cufftDoubleComplex *codeCorrNoFlipRaw_d;
        cufftDoubleComplex *codeCorrFlippedRaw_d;
        
        // The code-correlated signal reduced to one cp in size
        // Format: Prns * (S / N = number of samples per cp)
        cufftDoubleComplex *codeCorrNoFlipPerCP_d;
        cufftDoubleComplex *codeCorrFlippedPerCP_d;
        
        // Whether no-flip or flip gives the higher correlation score
        // Format: Prns 
        bool *noFlipIsLarger_d;
        
        
        // An array to hold the index of the next nav bit relative to this sample set
        // Format: one integer per channel (up to PRN_MAX)
        int *idxNextNavBit_d;
        
        // An array to hold the number of code periods elpased within this sample set
        // Format: one integer per channel (up to PRN_MAX)
        int *cpCompleted_d;
        
        // Ephemeris info
		std::vector<dsp::utils::ephSet_t> ephSetVec;
		int ephSetSize;
		dsp::utils::ephSet_t *ephSetArr;
		dsp::utils::ephSet_t *ephSetPtr_d;




        // Both intermediate and output
        dsp::utils::state_t<double> *currSatStates_d;




        // Outputs

        // Arrays holding the final output
        cufftDoubleComplex *codeCorrOut_d;
        cufftDoubleComplex *carrfftOut_d;
        cufftDoubleComplex *carrfftOutTrunc_d;
        double *txTime_d;
        bool *validPRN_d;

        // Debug: copies of the final output
        cufftDoubleComplex *codeCorrOut;
		cufftDoubleComplex *carrfftOut;

};

#endif

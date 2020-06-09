
#ifndef INC__DPINIT_H_
#define INC__DPINIT_H_


#include <stdint.h>
#include <cmath>
#include "dsp.h"
#include "module.h"
#include "auxil.h"
#include "ephhelper.h"
#include "consthelper.h"
#include "eigenwrapper.h"
#include "errorhandler.h"
//#include "constant.h"
#include "consthelper.h"
#include "rinexparse.h"
#include "statehelper.h"

#define DPI_STATE_SIZE 8

/** \brief DPInit module. */
class dsp::DPInit : public dsp::Module {

    public:

        DPInit();

        ~DPInit();

        /** \brief Any other initialization before running. */
        int Start(void* cuFlowStream);

        /** \brief Any de-initialization before de-constructor. */
        int Stop(void);
    
        /** \brief Update function to run every flow loop */
        int Update(void* cuFlowStream);     
        
        
        // Data parsing from file
        // Note: when editing this list, also edit the corresponding declaration in dpinit.cpp!!!!
        enum ScalarHandoffParams
        {
            rxTime,
            rxTime_a,
            X_ECEF,
            bytes_read,
            prn_list,
            rc,
            ri,
            fc,
            fi,
            cp,
            cp_timestamp,
            TOW,
            T_GD,
            C_uc,
            t_oe,
            t_oc,
            total,
            complete
        };   

    protected:
    
        /** \brief Get and parse fields from a file 
         * Adapted from https://stackoverflow.com/questions/12911299/read-csv-file-in-c*/
        void ParseField(char* line);
        
        double CorrectWeekCrossover(double time);
        int GetSatPos(dsp::utils::state_t<double> *satState, dsp::utils::eph_t *satEph, double txTime);
        void PerturbInitialization();

        /** Flag to signify if module is started. */
        char Started = 0;
    
        // Data to pass between modules
        
        double initrxTime;
        double initrxTime_a;

        long long int initByte;                          // To DPControl
        
        int initX_ECEF_allocated = DPI_STATE_SIZE;
        double* initX_ECEF = new double[DPI_STATE_SIZE];	// To EKF

        char *initPRN_d;                          // To DPControl
        char initPRNArr[CONST_PRN_MAX] = {(char)(-1)};
        int numPRNs;
        
        double *initRC_d;                          // To DPControl
        double initRCArr[CONST_PRN_MAX];
        
        double *initRI_d;                          // To DPControl
        double initRIArr[CONST_PRN_MAX];
        
        double *initFC_d;                          // To DPControl
        double initFCArr[CONST_PRN_MAX];
        
        double *initFI_d;                          // To DPControl
        double initFIArr[CONST_PRN_MAX];
        
        int *initCP_d;                          	// To DPControl
        int initCPArr[CONST_PRN_MAX];
        
        int *initCPTimestamp_d;                 	// To DPControl
        int initCPTimestampArr[CONST_PRN_MAX];
        
        std::vector<dsp::utils::ephSet_t> initEph;	// To SatPos

        int initTOE;								// To SatPos

        int initCPRefTOWArr[CONST_PRN_MAX];			// To SatPos

        double *initPArr;							// To cuEKF

        int initK;									// To cuEKF
        
        
        // Data parsing from file
        static const unsigned int HandoffFilenameCapacity = 1024;
        char HandoffFilename[HandoffFilenameCapacity] = {0};
        FILE* handoffParamsFile;
        auxil::EnumParser<ScalarHandoffParams> handoffParser;

        static const unsigned int RINEXFilenameCapacity = 1024;
        char RINEXFilename[RINEXFilenameCapacity] = {0};
        FILE* RINEXParamsFile;


        // Shifts in the initialization settings to ensure the grid isn't just sitting on the center point
        float initDeltaX;
        float initDeltaY;
        float initDeltaZ;
        float initDeltaT;
        double rangeDiff;


        // Debug
        int loopCounter = 0;


    private:
        /** \brief Flag to keep the thread loop running. */
        volatile bool KeepRunning = 0;

};

#endif

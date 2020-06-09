
#include <cuda_runtime.h>
#include <cuComplex.h>
#include <cufft.h>
#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
#include "dpinit.h"


/**
 * \brief The parameters and keys needed for DPE initialization
 *
 * Note: when editing this list, also edit the corresponding declaration in dpinit.h!!!!!
 */
namespace auxil {
template<>
EnumParser<dsp::DPInit::ScalarHandoffParams>::EnumParser() {
    enumMap["rxTime"]=dsp::DPInit::ScalarHandoffParams::rxTime;
    enumMap["rxTime_a"]=dsp::DPInit::ScalarHandoffParams::rxTime_a;
    enumMap["X_ECEF"]=dsp::DPInit::ScalarHandoffParams::X_ECEF;
    enumMap["bytes_read"]=dsp::DPInit::ScalarHandoffParams::bytes_read;
    enumMap["prn_list"]=dsp::DPInit::ScalarHandoffParams::prn_list;
    enumMap["rc"]=dsp::DPInit::ScalarHandoffParams::rc;
    enumMap["ri"]=dsp::DPInit::ScalarHandoffParams::ri;
    enumMap["fc"]=dsp::DPInit::ScalarHandoffParams::fc;
    enumMap["fi"]=dsp::DPInit::ScalarHandoffParams::fi;
    enumMap["cp"]=dsp::DPInit::ScalarHandoffParams::cp;
    enumMap["cp_timestamp"]=dsp::DPInit::ScalarHandoffParams::cp_timestamp;
    enumMap["TOW"]=dsp::DPInit::ScalarHandoffParams::TOW;
    enumMap["T_GD"]=dsp::DPInit::ScalarHandoffParams::T_GD;
    enumMap["C_uc"]=dsp::DPInit::ScalarHandoffParams::C_uc;
    enumMap["t_oe"]=dsp::DPInit::ScalarHandoffParams::t_oe;
    enumMap["t_oc"]=dsp::DPInit::ScalarHandoffParams::t_oc;
    enumMap["total"]=dsp::DPInit::ScalarHandoffParams::total;
    enumMap["complete"]=dsp::DPInit::ScalarHandoffParams::complete;
}
}



/**
 *   \brief Make adjustments to the initialized values as provided by parameters
 *
 *   Objective: if you want to make any adjustments to initialization values -- such as random starting states --
 *   adjust the parameters here.
 *
 */
void dsp::DPInit::PerturbInitialization() {

	initX_ECEF[0] += initDeltaX;
	initX_ECEF[1] += initDeltaY;
	initX_ECEF[2] += initDeltaZ;
	initX_ECEF[3] += initDeltaT;

}





dsp::DPInit::DPInit(){
    ModuleName = "DPInit";
    AllocateInputs(0); 	// DPInit bootstraps everything
    AllocateOutputs(14);

    ConfigOutput(0, "StartByte", INT_t, VALUE, HOST, 1, (void*)&initByte, 0);                                          // To DPControl

    ConfigOutput(1, "InitPRN", CHAR_t, VALUE, HOST, VECTORLENGTH_ANY, (void*)initPRNArr, 0);                            // To DPControl
    ConfigOutput(2, "InitCodePhase", DOUBLE_t, VALUE, HOST, VECTORLENGTH_ANY, (void*)initRCArr, 0);                     // To DPControl
    ConfigOutput(3, "InitCarrierPhase", DOUBLE_t, VALUE, HOST, VECTORLENGTH_ANY, (void*)initRIArr, 0);                  // To DPControl
    ConfigOutput(4, "InitCodeFrequency", DOUBLE_t, FREQUENCY_HZ, HOST, VECTORLENGTH_ANY, (void*)initFCArr, 0);          // To DPControl
    ConfigOutput(5, "InitCarrierFrequency", DOUBLE_t, FREQUENCY_HZ, HOST, VECTORLENGTH_ANY, (void*)initFIArr, 0);       // To DPControl
    ConfigOutput(6, "InitElapsedCodePeriods", INT_t, VALUE, HOST, VECTORLENGTH_ANY, (void*)initCPArr, 0);               // To DPControl
    ConfigOutput(7, "InitReferenceCodePeriods", INT_t, VALUE, HOST, VECTORLENGTH_ANY, (void*)initCPTimestampArr, 0);    // To DPControl
    ConfigOutput(8, "InitCPRefTOW", INT_t, VALUE, HOST, 1, (void*)&initCPRefTOWArr, 0);								// To BatchCorr

    ConfigOutput(9, "InitX", DOUBLE_t, STATE, HOST, VECTORLENGTH_ANY, (void*)&initX_ECEF, 0);                           // To EKF
    ConfigOutput(10, "InitP", DOUBLE_t, COVARIANCE, HOST, VECTORLENGTH_ANY, (void*)&initPArr, 0);                           // To EKF
    ConfigOutput(11, "InitK", INT_t, VALUE, HOST, 1, (void*)&initK, 0);                           // To EKF

    ConfigOutput(12, "InitRXTime", DOUBLE_t, VALUE, HOST, 1, (void*)&initrxTime, 0);								// To SatPos, DPMeas
    ConfigOutput(13, "InitEph", UNDEFINED_t, EPHEMS, HOST, 1, (void*)&initEph, 0);									// To SatPos

    
    std::clog << "[" << ModuleName << "] Configured outputs" << std::endl;
    
    InsertParam("HandoffFilename", (void*)&HandoffFilename, CHAR_t, HandoffFilenameCapacity, 0);
    InsertParam("RINEXFilename", (void*)&RINEXFilename, CHAR_t, RINEXFilenameCapacity, 0);
    
    InsertParam("InitDeltaX", (void*)&initDeltaX, FLOAT_t, sizeof(float), 1);
    InsertParam("InitDeltaY", (void*)&initDeltaY, FLOAT_t, sizeof(float), 1);
    InsertParam("InitDeltaZ", (void*)&initDeltaZ, FLOAT_t, sizeof(float), 1);
    InsertParam("InitDeltaT", (void*)&initDeltaT, FLOAT_t, sizeof(float), 1);

}


dsp::DPInit::~DPInit(){
    if (Started) Stop();
    delete [] expectedInputs;
    delete [] inputs;
    delete [] outputs;
}


/**
 *  Note: Don't do any cudaMalloc's in Start! This is not a .cu file!
 *  All outputs should be on HOST memory, and the modules are responsible
 *  for copying them to device for their use where needed.
 */
int dsp::DPInit::Start(void* cuFlowStream){
    if (Started) {
        std::clog << "[" << ModuleName << "] Start: Already Started." << std::endl;
        return 0;
    }
    std::clog << "[" << ModuleName << "] Starting ... " << std::flush;

    
    /**
     *  Parse Ephemerides
     */

    // Grab the file from the parameters
    RINEXParamsFile = fopen(this->RINEXFilename, "r");
    if (RINEXParamsFile == NULL) {
    	std::clog << "[" << ModuleName << "]" << " Open RINEXParamsFile failed: " << this->RINEXFilename << std::endl;
    	return -1;
    }

	// Rinex satellite system chooser
	//  (included for posterity compatibility with RTKLIB if multi-constellation support is desired in the future)
	const char *opt = "-SYS=G";
	if(dsp::utils::ReadRinex(RINEXParamsFile, opt, initEph) != 0) {
		std::clog << "[" << ModuleName << "] Stop: Failed to read in RINEX file" << std::endl;
		return -1;
	}


    /**
     *  Parse handoff params
     */
    handoffParamsFile = fopen(this->HandoffFilename, "r");
    if (handoffParamsFile == NULL) {
    	std::clog << "[" << ModuleName << "]" << " Open handoffParamsFile failed: " << this->HandoffFilename << std::endl;
    	return -1;
    }
    char curLine[1024];
    while (fgets(curLine, 1024, handoffParamsFile)) {
        char* tmp = strdup(curLine);
        ParseField(tmp); // NOTE: strtok clobbers tmp
        free(tmp);
    }

    // Also initialize the P matrix (until there's a better way to set this up)
    initPArr = auxil::MakeIMatrix<double>(DPI_STATE_SIZE, DPI_STATE_SIZE);
    // And start at the first measurement
    initK = 0;



    /**
     *  Adjust initialization
     */
    PerturbInitialization();


    // Update the VectorLengths now that that's known
    errCheckModSt(UpdateOutput(1, 	numPRNs, 	(void*)initPRNArr, 0));
    errCheckModSt(UpdateOutput(2, 	numPRNs, 	(void*)initRCArr, 0));
    errCheckModSt(UpdateOutput(3, 	numPRNs, 	(void*)initRIArr, 0));
    errCheckModSt(UpdateOutput(4, 	numPRNs, 	(void*)initFCArr, 0));
    errCheckModSt(UpdateOutput(5, 	numPRNs, 	(void*)initFIArr, 0));
    errCheckModSt(UpdateOutput(6, 	numPRNs, 	(void*)initCPArr, 0));
    errCheckModSt(UpdateOutput(7, 	numPRNs, 	(void*)initCPTimestampArr, 0));
    errCheckModSt(UpdateOutput(8, 	numPRNs, 	(void*)&initCPRefTOWArr, 0));

    errCheckModSt(UpdateOutput(9, 	initX_ECEF_allocated, (void*)initX_ECEF, 0));
    errCheckModSt(UpdateOutput(10, 	initX_ECEF_allocated*initX_ECEF_allocated, (void*)initPArr, 0));
    errCheckModSt(UpdateOutput(11, 	1, (void*)&initK, 0));

    errCheckModSt(UpdateOutput(12, 	1, (void*)&initrxTime, 0));
    errCheckModSt(UpdateOutput(13, 	initEph.size(), (void*)&initEph, 0));



    std::clog << "[" << ModuleName << "] Updated outputs" << std::endl;



    Started = 1;

    std::clog << "Started." << std::endl;
    return 0;
}

int dsp::DPInit::Stop(){
    int ret = 0;

    if (Started == 0) {
        std::clog << "[" << ModuleName << "] Stop: Wasn't Started." << std::endl;
        return 0;
    }
    std::clog << "[" << ModuleName << "] Stopping ... " << std::flush;
    
    
    

    Started = 0;
    std::clog << "Stopped." << std::endl;

    return ret;
}


int dsp::DPInit::Update(void* cuFlowStream){

	if (loopCounter % 500 == 0) {
		std::clog << "[" << ModuleName << "] Started iteration " << loopCounter << std::endl;
	}
    loopCounter++;


    // Bootleg way to end the pipeline after a certain number of iterations
    if (loopCounter >= 3000) {
    	return -1;
    	//Stop();
    	//KeepRunning = false;
    }

    return 0;
}


/**
 * \brief Grab the parameters from the current line from the initialization file
 *
 * \param line	The line of initialization parameters to be processed
 *
 */
void dsp::DPInit::ParseField(char* line) {
    const char* fieldType;
    const char* tok;
    
    int i = 0;
    
    fieldType = strtok(line, ",");
    
    
    switch (handoffParser.ParseSomeEnum(fieldType)) {

		// GPS time at initialization
		case rxTime:
			tok = strtok(NULL, "\r");
			initrxTime = atof(tok);
			//std::cout << "[DPInit] Found rxTime: " << initrxTime << std::endl;
			break;

		// GPS adjusted time at initialization
		case rxTime_a:
			tok = strtok(NULL, "\r");
			initrxTime_a = atof(tok);
			//std::cout << "[DPInit] Found rxTime_a: " << initrxTime_a << std::endl;
			break;

		// Initial state
		case X_ECEF:
			//std::cout << "[DPInit] Found X_ECEF: ";
			for(tok = strtok(NULL, ","); tok && *tok; tok = strtok(NULL, ",")) {
				if(initX_ECEF_allocated < i) {
					if (!(initX_ECEF = (double *)realloc(initX_ECEF,sizeof(double)*(++initX_ECEF_allocated)))) {
						std::cout << "[" << ModuleName << "] X_ECEF initialization malloc error: " << std::endl;
						free(initX_ECEF);
						return;
					}
				}
				double curX_ECEF = atof(tok);
				initX_ECEF[i] = curX_ECEF;
				//std::cout << initX_ECEFArr[i] << " ";
				i++;
			}
			break;

    	// File read info
        case bytes_read:
            tok = strtok(NULL, "\r");
            initByte = atoll(tok); // Need long long int to be able to read past 214 secs for a 2500kHz file
            //std::cout << "[DPInit] Found bytes_read: " << startByte << std::endl;
            break;


        // Locked satellites info
        case prn_list:
        	//std::cout << "[DPInit] Found prn_list: ";
            for(tok = strtok(NULL, ","); tok && *tok; tok = strtok(NULL, ",")) {
                char curPRN = (char)atoi(tok);
            	initPRNArr[i] = curPRN;
                //std::cout << static_cast<int>(initPRNArr[i]) << " ";
                i++;
            }
            //std::cout << std::endl;
            numPRNs = i;
            break;
        
        // Channel params
        case rc:
        	//std::cout << "[DPInit] Found rc: ";
            for(tok = strtok(NULL, ","); tok && *tok; tok = strtok(NULL, ",")) {
                double curRC = atof(tok);
            	initRCArr[i] = curRC;
                //std::cout << initRCArr[i] << " ";
                i++;
            }
            //std::cout << std::endl;
            break;
        case ri:
        	//std::cout << "[DPInit] Found ri: ";
            for(tok = strtok(NULL, ","); tok && *tok; tok = strtok(NULL, ",")) {
                double curRI = atof(tok);
            	initRIArr[i] = curRI;
                //std::cout << initRIArr[i] << " ";
                i++;
            }
            //std::cout << std::endl;
            break;
        case fc:
        	//std::cout << "[DPInit] Found fc: ";
            for(tok = strtok(NULL, ","); tok && *tok; tok = strtok(NULL, ",")) {
                double curFC = atof(tok);
            	initFCArr[i] = curFC;
                //std::cout << initFCArr[i] << " ";
                i++;
            }
            //std::cout << std::endl;
            break;
        case fi:
        	//std::cout << "[DPInit] Found fi: ";
            for(tok = strtok(NULL, ","); tok && *tok; tok = strtok(NULL, ",")) {
                double curRI = atof(tok);
            	initFIArr[i] = curRI;
                //std::cout << initFIArr[i] << " ";
                i++;
            }
            //std::cout << std::endl;
            break;
        case cp:
        	//std::cout << "[DPInit] Found cp: ";
            for(tok = strtok(NULL, ","); tok && *tok; tok = strtok(NULL, ",")) {
                int curCP = atof(tok);
            	initCPArr[i] = curCP;
                //std::cout << initCPArr[i] << " ";
                i++;
            }
            //std::cout << std::endl;
            break;
        case cp_timestamp:
        	//std::cout << "[DPInit] Found cp_timestamp: ";
            for(tok = strtok(NULL, ","); tok && *tok; tok = strtok(NULL, ",")) {
                int curCPTimestamp = atof(tok);
            	initCPTimestampArr[i] = curCPTimestamp;
                //std::cout << initCPTimestampArr[i] << " ";
                i++;
            }
            //std::cout << std::endl;
            break;

        // Grab the Time of Ephemeris used by the DP initialization
        case t_oe:
        	tok = strtok(NULL, ",");
        	initTOE = atoi(tok);
        	break;

		// Grab the TOW with which all code periods are referenced
		case TOW:
        	//std::cout << "[DPInit] Found TOW: ";
            for(tok = strtok(NULL, ","); tok && *tok; tok = strtok(NULL, ",")) {
                int curCPRefTOW = atof(tok);
            	initCPRefTOWArr[i] = curCPRefTOW;
                //std::cout << initTOWArr[i] << " ";
                i++;
            }
            //std::cout << std::endl;
			break;

        
        // Default
        default:
            for(tok = strtok(NULL, ","); tok && *tok; tok = strtok(NULL, ",")) {
                // Don't know what to do with this field, so advance to the end of it
            }
            //std::cout << "[DPInit] Found unknown, advancing through this line." << std::endl;
            break;
	}
}

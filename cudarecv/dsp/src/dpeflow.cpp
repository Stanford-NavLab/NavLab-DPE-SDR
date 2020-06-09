
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <fcntl.h>
#include <unistd.h>
#include <sched.h>
#include <pthread.h>
#include <ctime>
#include <stdlib.h>
#include <string.h>
#include "dpeflow.h"
#include "dpinit.h"
#include "sampleblock.h"
#include "batchcorrscores.h"
#include "batchcorrmanifold.h"
#include "cuekf.h"
#include "cuchanmgr.h"
#include "datalogger.h"
#include "satpos.h"
#include "gridhelper.h"



int dsp::DPEFlow::LoadFlow(const char *filename) {

	// Files for the algorithm
	char* homePath = getenv("HOME");

    char rawFilename[] 	= "/Desktop/demofiles/static_opensky_20180705_190000_usrp6_2500kHz.dat";
    char dummyRawFilename[256];
    strcpy(dummyRawFilename, homePath);
    strcat(dummyRawFilename, rawFilename);

    char handoffFilename[] = "/Desktop/demofiles/handoff_params_usrp6.csv";
    char dummyHandoffFilename[256];
    strcpy(dummyHandoffFilename, homePath);
    strcat(dummyHandoffFilename, handoffFilename);

    char RINEXFilename[]	= "/Desktop/demofiles/nist1860.18n";	// Consider automating this as CDDIS ftp download
    char dummyRINEXFilename[256];
    strcpy(dummyRINEXFilename, homePath);
    strcat(dummyRINEXFilename, RINEXFilename);

    char gridFilename[] 	= "/Desktop/demofiles/rngrid3.csv";
    char dummyGridFilename[256];
    strcpy(dummyGridFilename, homePath);
    strcat(dummyGridFilename, gridFilename);




    // Configure the flow
    Mods.reserve(7);    // Number of Modules in this flow
    Mods.push_back(new dsp::DPInit);
    Mods.push_back(new dsp::SampleBlock);
    Mods.push_back(new dsp::BatchCorrScores);
    Mods.push_back(new dsp::BatchCorrManifold);
    Mods.push_back(new dsp::cuEKF);
    Mods.push_back(new dsp::cuChanMgr);
    Mods.push_back(new dsp::DataLogger("XECEFLogger"));



    // Set parameters to different modules
    errCheck(SetModParam("SampleBlock", "SamplingFrequency", 2.5e6)); // TODO: have this derived by file name!
    errCheck(SetModParam("SampleBlock", "RunLive", false));
    errCheck(SetModParam("SampleBlock", "Filename", dummyRawFilename));
    errCheck(SetModParam("DPInit", "HandoffFilename", dummyHandoffFilename));
    errCheck(SetModParam("DPInit", "RINEXFilename", dummyRINEXFilename));

    // Used to push the initialization of DPInit around to make sure that the grid doesn't just sit on the first point
    errCheck(SetModParam("DPInit", "InitDeltaX", (float)0.0));
    errCheck(SetModParam("DPInit", "InitDeltaY", (float)0.0));
    errCheck(SetModParam("DPInit", "InitDeltaZ", (float)0.0));
    errCheck(SetModParam("DPInit", "InitDeltaT", (float)0.0));

    // Does SampleBlock need T = 0.001?
    double T = 0.02; // Snippet of samples size
    errCheck(SetModParam("SampleBlock", "SampleLength", T)); // THINGS WILL BREAK IF T IS CHANGED!!!!!!
    errCheck(SetModParam("cuEKF", "SampleLength", T)); // THINGS WILL BREAK IF T IS CHANGED!!!!!!
    errCheck(SetModParam("BatchCorrManifold", "PosGridDimSize", (int)25));
    errCheck(SetModParam("BatchCorrManifold", "VelGridDimSize", (int)25));
    errCheck(SetModParam("BatchCorrManifold", "GridDimSpacing", (float)1.0)); // This can be overridden in BatchCorrManifold
    errCheck(SetModParam("BatchCorrManifold", "GridType", dsp::utils::ManifoldGridTypes::Uniform));
    errCheck(SetModParam("BatchCorrManifold", "LPower", (int)1));
    errCheck(SetModParam("cuChanMgr", "DopplerSign", (int)1));	// Either +1 or -1
    
    errCheck(SetModParam("cuEKF", "EnableEKF", (bool)false));


    // Name the output file from the current date and time
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];

    time (&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(buffer,sizeof(buffer),"%d-%m-%Y_%H-%M-%S",timeinfo);
    std::string str(buffer);
    // Helpers for the name of the output file
    std::string XFileNamePrefix = "/home/ubuntu/output/inves/";
    std::string XFileNameSuffix = "-XFile.csv";
    std::string GridNameSuffix = "-Grid.csv";
    std::string rcSuffix = "-rc.csv";
    std::string fcSuffix = "-fc.csv";
    std::string riSuffix = "-ri.csv";
    std::string fiSuffix = "-fi.csv";
    std::string cpSuffix = "-cp.csv";

    std::string XFileName;
    XFileName = XFileNamePrefix + str + XFileNameSuffix;
    std::string rcFileName;
	rcFileName = XFileNamePrefix + str + rcSuffix;
    std::string fcFileName;
	fcFileName = XFileNamePrefix + str + fcSuffix;
    std::string riFileName;
	riFileName = XFileNamePrefix + str + riSuffix;
    std::string fiFileName;
	fiFileName = XFileNamePrefix + str + fiSuffix;
    std::string cpFileName;
	cpFileName = XFileNamePrefix + str + cpSuffix;
	std::string GridFileName;
    GridFileName = XFileNamePrefix + str + GridNameSuffix;

    errCheck(SetModParam("XECEFLogger", "Filename", XFileName.c_str()));
    errCheck(SetModParam("XECEFLogger", "CSV", true));


    errCheck(SetModParam("BatchCorrManifold", "GridLogFileName", GridFileName.c_str()));
    errCheck(SetModParam("BatchCorrManifold", "LoadPosGridFilename", dummyGridFilename));



    // Connect Ports -- grouped by source module
    
    // DPInit gives the starting sample and starting channel params
    errCheck(ConnectPort("DPInit", "StartByte", 				"SampleBlock", 			"StartByte"));
    errCheck(ConnectPort("DPInit", "InitX", 					"cuEKF", 				"InitX"));
    errCheck(ConnectPort("DPInit", "InitP", 					"cuEKF", 				"InitP"));
    errCheck(ConnectPort("DPInit", "InitK", 					"cuEKF", 				"InitK"));
    errCheck(ConnectPort("DPInit", "InitEph", 					"cuChanMgr", 			"InitEph"));
    errCheck(ConnectPort("DPInit", "InitPRN", 					"cuChanMgr", 			"InitPRN"));
    errCheck(ConnectPort("DPInit", "InitCodePhase", 			"cuChanMgr", 			"InitCodePhase"));
    errCheck(ConnectPort("DPInit", "InitCarrierPhase", 			"cuChanMgr", 			"InitCarrierPhase"));
    errCheck(ConnectPort("DPInit", "InitCodeFrequency", 		"cuChanMgr", 			"InitCodeFrequency"));
    errCheck(ConnectPort("DPInit", "InitCarrierFrequency", 		"cuChanMgr", 			"InitCarrierFrequency"));
    errCheck(ConnectPort("DPInit", "InitElapsedCodePeriods", 	"cuChanMgr", 			"InitElapsedCodePeriods"));
    errCheck(ConnectPort("DPInit", "InitReferenceCodePeriods", 	"cuChanMgr", 			"InitReferenceCodePeriods"));
    errCheck(ConnectPort("DPInit", "InitCPRefTOW", 				"cuChanMgr", 			"InitCPRefTOW"));
    errCheck(ConnectPort("DPInit", "InitRXTime", 				"cuChanMgr", 			"InitRXTime"));

    // SampleBlock outputs samples
    errCheck(ConnectPort("SampleBlock", "Samples", 				"BatchCorrScores", 		"Samples"));
    errCheck(ConnectPort("SampleBlock", "SamplingFrequency", 	"BatchCorrScores", 		"SamplingFrequency"));
    errCheck(ConnectPort("SampleBlock", "SampleLength",			"BatchCorrScores", 		"SampleLength"));
    errCheck(ConnectPort("SampleBlock", "SamplingFrequency", 	"BatchCorrManifold", 	"SamplingFrequency"));
    errCheck(ConnectPort("SampleBlock", "SampleLength", 		"BatchCorrManifold", 	"SampleLength"));
    errCheck(ConnectPort("SampleBlock", "SampleLength",			"cuChanMgr", 			"SampleLength"));
    
    // BatchCorrScores gives a vector of code and carrier scores
    errCheck(ConnectPort("BatchCorrScores", "CodeScores", 		"BatchCorrManifold", "CodeScores"));
    errCheck(ConnectPort("BatchCorrScores", "CarrScores", 		"BatchCorrManifold", "CarrScores"));
    errCheck(ConnectPort("BatchCorrScores", "NumFFTPoints",		"BatchCorrManifold", "NumFFTPoints"));

    // ChannelManager keeps track of the channel parameters and calculates the corresponding satellite positions
    errCheck(ConnectPort("cuChanMgr", "CodePhaseStart", 	"BatchCorrScores", "CodePhaseStart"));
    errCheck(ConnectPort("cuChanMgr", "CodeFrequency", 		"BatchCorrScores", "CodeFrequency"));
    errCheck(ConnectPort("cuChanMgr", "CarrierPhaseStart", 	"BatchCorrScores", "CarrierPhaseStart"));
	errCheck(ConnectPort("cuChanMgr", "CarrierFrequency", 	"BatchCorrScores", "CarrierFrequency"));
	errCheck(ConnectPort("cuChanMgr", "cpReference", 		"BatchCorrScores", "cpReference"));
	errCheck(ConnectPort("cuChanMgr", "cpElapsedStart", 	"BatchCorrScores", "cpElapsedStart"));
	errCheck(ConnectPort("cuChanMgr", "DopplerSign", 		"BatchCorrScores", "DopplerSign"));
	errCheck(ConnectPort("cuChanMgr", "ValidPRNs", 			"BatchCorrScores", "ValidPRNs"));

    errCheck(ConnectPort("cuChanMgr", "CodeFrequency", 		"BatchCorrManifold", "CodeFrequency"));
    errCheck(ConnectPort("cuChanMgr", "CarrierFrequency", 	"BatchCorrManifold", "CarrierFrequency"));
    errCheck(ConnectPort("cuChanMgr", "rxTime",				"BatchCorrManifold", "rxTime"));
    errCheck(ConnectPort("cuChanMgr", "txTime", 			"BatchCorrManifold", "txTime"));
    errCheck(ConnectPort("cuChanMgr", "DopplerSign",		"BatchCorrManifold", "DopplerSign"));
    errCheck(ConnectPort("cuChanMgr", "SatStates",			"BatchCorrManifold", "SatStates"));
    errCheck(ConnectPort("cuChanMgr", "ENU2ECEFMat", "BatchCorrManifold", "ENU2ECEFMat"));
    errCheck(ConnectPort("cuChanMgr", "SatStatesOld", "BatchCorrManifold", "SatStatesOld"));
    
    errCheck(ConnectPort("cuChanMgr", "CodePhaseEnd", 		"BatchCorrManifold", "CodePhase"));
	errCheck(ConnectPort("cuChanMgr", "CarrierPhaseEnd", 	"BatchCorrManifold", "CarrierPhase"));
	errCheck(ConnectPort("cuChanMgr", "cpRefTOW", 			"BatchCorrManifold", "cpRefTOW"));
	errCheck(ConnectPort("cuChanMgr", "cpRef", 				"BatchCorrManifold", "cpRef"));
	errCheck(ConnectPort("cuChanMgr", "cpElapsedEnd", 		"BatchCorrManifold", "cpElapsedEnd"));

	// Debug -- if you want to view channel parameters through the processing
	//errCheck(ConnectPort("cuChanMgr", "CodePhaseStart", 	"rcLogger", "Data"));
	//errCheck(ConnectPort("cuChanMgr", "CodeFrequency", 		"fcLogger", "Data"));
	//errCheck(ConnectPort("cuChanMgr", "CarrierPhaseStart", 	"riLogger", "Data"));
	//errCheck(ConnectPort("cuChanMgr", "CarrierFrequency", 	"fiLogger", "Data"));
	//errCheck(ConnectPort("cuChanMgr", "cpElapsedStart", 	"cpLogger", "Data"));

    // BatchCorrManifold gives a measurement to the Kalman filter
    errCheck(ConnectPort("BatchCorrManifold", "zVal", "cuEKF", "zVal"));
    errCheck(ConnectPort("BatchCorrManifold", "RVal", "cuEKF", "RVal"));
    errCheck(ConnectPort("BatchCorrManifold", "TimeGrid", "cuChanMgr", "TimeGrid"));


    // The Kalman filter returns the current estimate
    //     Manifold should be centered on the best estimate of the EKF, taking the time estimate
    //     Scores could use the time estimate also, but, to match PyGNSS, will take the measurement estimate
    //         since it should be safe to assume the code/carr frequencies won't change much between updates
    errCheck(ConnectPort("cuEKF", "xCurrk1k1", "cuChanMgr", 		"xCurrk1k1"));
    errCheck(ConnectPort("cuEKF", "xCurrkk1",  "cuChanMgr", 		"xCurrkk1"));
    errCheck(ConnectPort("cuEKF", "xCurrkk1",  "BatchCorrManifold", "xCurrkk1"));
    errCheck(ConnectPort("cuEKF", "xCurrk1k1", "XECEFLogger", 	"Data"));


    
    
    // Report success to the user
    std::clog << "[DPEFlow] Completed LoadFlow." << std::endl;

    return 0;
}



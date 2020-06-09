#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <pthread.h>
#include <signal.h>
#include "startoptions.h"
#include "console.h"
#include "cmdParser.h"
#include "flowmgr.h"

#include "dsp.h"
//#include "acquisition.h"
#include "dpeflow.h"

// For math functions used by initialization shift
#include <math.h>
#include <stdlib.h>

using namespace std;
console::CmdParser * cmdMgr = NULL;
static dsp::FlowMgr* fm = nullptr;

//void sigHandler(int sig);
volatile char KeepRunning = 1;

void sigHandler(int sig){
    signal(sig, SIG_IGN);
    /*if (fm && fm->EmergencyStop()){
        std::cout << "[SIGINT] Press Ctrl+C again to quit the program." << std::endl;
    } else {*/
        cout << "[SIGINT] Signal Caught... exiting" << endl;
        KeepRunning = 0;
    //}
}

/** \brief  Main entry to program. */
int main(int argc, char** argv) {

    std::cout << "GNSS Receiver by Gao Research Group" <<endl;
    std::cout << "Version 0.0 Alpha" <<endl;
    //dsp::testAcquisition();
    StartOptions options(argc, argv);

	if (options.error) return -1;

    //cout << "console: " << options.console << endl;
    //if (options.fromFile)
    //    cout << "from file: " << options.filename << endl;

    // Setup POSIX signals to handle ctrl+C ourselves
    KeepRunning = 1;

    struct sigaction sigIntHandler;
    sigIntHandler.sa_handler = sigHandler;
    sigemptyset(&sigIntHandler.sa_mask);
    sigIntHandler.sa_flags = 0;
    sigaction(SIGINT, &sigIntHandler, NULL);

    /* Set main thread as the lowest priority in
     * real time round robin policy. */
    pthread_t mainThread = pthread_self();
    struct sched_param param;
    param.sched_priority = sched_get_priority_min(SCHED_RR);
    pthread_setschedparam(mainThread, SCHED_RR, &param);




    // Primary way of running CUDARecv -- ask the user what to do
    // Follow the README for instructions on how to interact with the program
    if (options.console) {
        cmdMgr = new console::CmdParser("\033[01;32mgnss\033[00m$ ");
        fm = new dsp::FlowMgr;
        console::initCommonCmd(cmdMgr);
        console::initFlowCmd(cmdMgr,fm);

        while ((cmdMgr->execOneCmd()) && KeepRunning);

        delete fm;
        delete cmdMgr;
    } else

	return 0;





    // Automated method 2: load and start -- don't bother interacting with the user, just launch DPEFlow
    /*
	dsp::DPEFlow myFlow;
	int i = myFlow.LoadFlow(NULL);
	if (i == 0) myFlow.Start();
	//while ((cmdMgr->execOneCmd()) && KeepRunning);
	while (!myFlow.CheckFlowState()) {
		sleep(1);
	}
	myFlow.Stop();
	cout << "[MAIN] Stopped flow." << endl;

    */



    // Automated method 2: random initial receiver state guess -- perturb the initial state and run repeatedly, indexing the output files
    /*
    // Rotate the offset into ECEF frame (note: this only works for the usrp 6 simulated datasets!)
    double r_enu2ecef[9] = {0.9995, -0.0199, 0.0237, 0.0309, 0.6440, -0.7644, 0, 0.7648, 0.6443}; // Receiver in Urbana-Champaign, IL, USA

    // Name the output file from the current date and time
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];


    // Set 1
    time (&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(buffer,sizeof(buffer),"%d-%m-%Y_%H-%M-%S",timeinfo);
    std::string fnameTime(buffer);
    // Helpers for the name of the output file
    std::string XFileNamePrefix = "/home/ubuntu/output/test/";
    std::string XFileNameSuffix = "-XFile.csv";
    std::string XFileName;


    cmdMgr = new console::CmdParser("\033[01;32mgnss\033[00m$ ");
    console::initCommonCmd(cmdMgr);
    console::initFlowCmd(cmdMgr,fm);


    srand(time(0));
	ofstream shiftFile;
	std::string shiftFileName;
	shiftFileName = XFileNamePrefix + "shiftFile_" + fnameTime + ".csv";
	shiftFile.open(shiftFileName);


    for (int idx = 0; idx < 100; idx++) {

		dsp::DPEFlow myFlow;
		int i = myFlow.LoadFlow(NULL);


		// Specify how much spread to use
		//float shiftMin = 50.0;
		//float shiftMax = 80.0;

		float shiftRange = 30.0;
		//float shiftRange = 25.0;
		float shiftBottom = 50.0;
		//float shiftBottom = 25.0;


		// (New) Horizontal shift
		float magShift = static_cast<float>(rand()) / (static_cast<float>(RAND_MAX/(2*shiftRange)))-shiftRange;
		if (magShift >= 0) {
			magShift += shiftBottom;
		}
		else {
			magShift -= shiftBottom;
		}
		float thetaShift = static_cast<float>(rand()) / (static_cast<float>(RAND_MAX/(3.1415926535*2)));
		float yShiftENU = magShift * sin(thetaShift);
		float xShiftENU = magShift * cos(thetaShift);


    	// Vertical shift
    	float zShiftENU = static_cast<float>(rand()) / (static_cast<float>(RAND_MAX/(2*shiftRange)))-shiftRange;

    	if (zShiftENU >= 0) {
			zShiftENU += shiftBottom;
		}
		else {
			zShiftENU -= shiftBottom;
		}


    	// Time shift
    	float tShift = static_cast<float>(rand()) / (static_cast<float>(RAND_MAX/(2*shiftRange)))-shiftRange;

    	if (tShift >= 0) {
			tShift += shiftBottom;
		}
		else {
			tShift -= shiftBottom;
		}

    	float xShift = r_enu2ecef[0] * xShiftENU + r_enu2ecef[1] * yShiftENU + r_enu2ecef[2] * zShiftENU;
    	float yShift = r_enu2ecef[3] * xShiftENU + r_enu2ecef[4] * yShiftENU + r_enu2ecef[5] * zShiftENU;
    	float zShift = r_enu2ecef[6] * xShiftENU + r_enu2ecef[7] * yShiftENU + r_enu2ecef[8] * zShiftENU;

    	//xShift = 0;
    	//yShift = 0;
    	//zShift = 0;
    	tShift = 0;

		myFlow.SetModParam("DPInit", "InitDeltaX", xShift);
		myFlow.SetModParam("DPInit", "InitDeltaY", yShift);
		myFlow.SetModParam("DPInit", "InitDeltaZ", zShift);
		myFlow.SetModParam("DPInit", "InitDeltaT", tShift);
		myFlow.SetModParam("BatchCorrManifold", "LoadPosGrid", (bool)false);

    	shiftFile << idx << ", " << xShift << ", " << yShift << ", " << zShift << ", " << tShift << std::endl;

		XFileName = XFileNamePrefix + fnameTime + "_idx" + std::to_string(idx) + XFileNameSuffix;

		myFlow.SetModParam("XECEFLogger", "Filename", XFileName.c_str());


		if (i == 0) myFlow.Start();
		//while ((cmdMgr->execOneCmd()) && KeepRunning);
		while (!myFlow.CheckFlowState()) {
			sleep(1);
		}
		myFlow.Stop();
		cout << "[MAIN] Stopped flow." << endl;
    }

    shiftFile.close();

     */




    // Automated method 3: resizable grids -- change grid spacing and run repeatedly, indexing the output files
    /*
    // Name the output file from the current date and time
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];

    time (&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(buffer,sizeof(buffer),"%d-%m-%Y_%H-%M-%S",timeinfo);
    std::string fnameTime(buffer);
    // Helpers for the name of the output file
    std::string XFileNamePrefix = "/home/ubuntu/output/bulk/";
    std::string XFileNameSuffix = "-XFile.csv";
    std::string XFileName;


    cmdMgr = new console::CmdParser("\033[01;32mgnss\033[00m$ ");
    console::initCommonCmd(cmdMgr);
    console::initFlowCmd(cmdMgr,fm);


    srand(static_cast<unsigned> (time(0)));
	ofstream shiftFile;
	std::string shiftFileName;
	shiftFileName = XFileNamePrefix + "shiftFile_" + fnameTime + ".csv";
	shiftFile.open(shiftFileName);

    for (int idx = 0; idx < 7; idx++) {

		dsp::DPEFlow myFlow;
		int i = myFlow.LoadFlow(NULL);

		XFileName = XFileNamePrefix + fnameTime + "_idx" + std::to_string(idx) + XFileNameSuffix;

		myFlow.SetModParam("XECEFLogger", "Filename", XFileName.c_str());

		myFlow.SetModParam("BatchCorrManifold", "GridDimSpacing", (float)((idx+1)*0.5 + 6.5));
		//myFlow.SetModParam("BatchCorrManifold", "GridDimSpacing", (float)((idx+1)+7));


		if (i == 0) myFlow.Start();
		//while ((cmdMgr->execOneCmd()) && KeepRunning);
		while (!myFlow.CheckFlowState()) {
			sleep(1);
		}
		myFlow.Stop();
		cout << "[MAIN] Stopped flow." << endl;
    }

    shiftFile.close();
     */



    return 0;

}

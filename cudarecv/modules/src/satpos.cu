
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
//#include "cuPrintf.cu"
#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <stdint.h>
#include "satpos.h"
#include "errorhandler.h"
#include "auxil.h"


// TODO: Get cuPrintf.cu working so that satellite position calculation error can be reported back

__device__ double
SPS_Correct_Week_Crossover(double time) {
	if (time > 302400.0) 		{ return time-604800.0;	}
	else if (time < -302400.00)	{ return time+604800.0;	}
	else						{ return time;			}
}



// Function to be called within Update_Sat_Pos to prevent function saturation
// Note: this is hardcoded to compute GPS ephemerides only!
// See ephemeris.c in RTKLIB for potentially useful code for multi-constellation
//
// NOTE: IF THERE ARE ERRORS, "FMOD" is defined differently than np.mod! Check those lines!!!
//
__device__ int
SPS_Get_Sat_Pos(dsp::utils::state_t<double> *satState, dsp::utils::eph_t *satEph, double txTime, int curBatch, int curSat, int curPos, int ephToe) {

	// Corrected mean motion
	double n = sqrt(dsp::utils::MU_GPS/(SATPOS_CUBE(satEph->A))) + satEph->deln;

	// Compute satellite clock corrections
	double tc = SPS_Correct_Week_Crossover(txTime - satEph->tocs); 					// Without corrections
	double clkb = satEph->f2*tc*tc + satEph->f1*tc + satEph->f0 - satEph->tgd[0];	// Without relativistic correction
	double tk = SPS_Correct_Week_Crossover(txTime - clkb - satEph->toes);			// Without relativistic correction

	// Mean anomaly
	double E, M, f, dfdE, dE = 1;
	E = M = fmod((satEph->M0 + n*tk), CONST_2PI);
	// Eccentric anomaly
	for (int eccIdx = 0; eccIdx < SATPOS_MAX_ITER_KEPLER && abs(dE) > SATPOS_RTOL_KEPLER; eccIdx++) {
		f = M - E + satEph->e * sin(E);
		dfdE = -1.0 + satEph->e * cos(E);
		dE = -f / dfdE;
		E = fmod(E + dE, CONST_2PI);
		//if (abs(dE) < SATPOS_RTOL_KEPLER) { break; } // Break here if convergence is achieved!
	}
	if (abs(dE) > SATPOS_RTOL_KEPLER) { return -1;	}

	// Add in relativistic corrections and compute clock drift
	double dtr = CONST_F*(satEph->e)*(satEph->sqrt_A)*sin(E);
	tc = txTime - (clkb + dtr) - satEph->tocs;
	clkb = satEph->f2*tc*tc + satEph->f1*tc + satEph->f0 + dtr - satEph->tgd[0];
	double clkd = satEph->f1 + 2.0*satEph->f2*tc;


	// Recompute tk with relativisitic correction
	tk = SPS_Correct_Week_Crossover(txTime - clkb - satEph->toes);

	// Recompute mean anomaly with relativisitic correction
	dE = 1;
	E = M = fmod((satEph->M0 + n*tk), (CONST_2PI));
	// Eccentric anomaly
	for (int eccIdx = 0; eccIdx < SATPOS_MAX_ITER_KEPLER && abs(dE) > SATPOS_RTOL_KEPLER; eccIdx++) {
		f = M - E + satEph->e * sin(E);
		dfdE = -1.0 + satEph->e * cos(E);
		dE = -f / dfdE;
		E = fmod(E + dE, CONST_2PI);
		//if (abs(dE) < SATPOS_RTOL_KEPLER) { break; } // Break here if convergence is achieved!
	}
	if (abs(dE) > SATPOS_RTOL_KEPLER) { return -1;	}

	// Compute helpers
	double sinE = sin(E);
	double cosE = cos(E);


	// True anomaly
	double v = atan2(sqrt(1.0-satEph->e_sqr) * sinE / (1.0-satEph->e*cosE), (cosE-satEph->e)/(1.0-satEph->e*cosE));
	// Argument of latitude
	double u = fmod(v + satEph->omg, CONST_2PI); // Don't need mod here? (Computing for guarantee)
	//double u = v + satEph->omg;

	// Second harmonic perturbations
	double cos2u = cos(2.0*u);
	double sin2u = sin(2.0*u);

	// Argument of latitude correction -> corrected argument of latitude
	u += satEph->cuc * cos2u + satEph->cus * sin2u;
	// Radius correction -> corrected radius
	double r = satEph->A * (1.0 - satEph->e * cosE)   +   satEph->crc * cos2u + satEph->crs * sin2u;
	// Orbital inclination correction -> corrected inclination
	double i = satEph->i0 + satEph->idot * tk   +   satEph->cic * cos2u + satEph->cis * sin2u;

	// Corrected longitude of node
	double omegak = fmod(satEph->OMG0 + (satEph->OMGd-CONST_OEDot)*tk - CONST_OEDot*satEph->toes, CONST_2PI);
	//double omegak = satEph->OMG0 + (satEph->OMGd - CONST_OEDot)*tk - CONST_OEDot * satEph->toes; // Don't need mod here?

	// Positions in orbital plane
	double x_op = r * cos(u);
	double y_op = r * sin(u);

	double cos_omegak = cos(omegak);
	double sin_omegak = sin(omegak);
	double cosi = cos(i);
	double sini = sin(i);

	// Assign position states
	double state_x = x_op * cos_omegak - y_op * sin_omegak * cosi;
	satState->x = state_x;
	double state_y = x_op * sin_omegak + y_op * cos_omegak * cosi;
	satState->y = state_y;
	double state_z = y_op * sini;
	satState->z = state_z;

	satState->delta_t = clkb;



	// Velocity calculation

	// Second harmonic perturbations
	cos2u = cos(2.0*u);
	sin2u = sin(2.0*u);

	double edot = n / (1.0 - satEph->e*cosE);
	double vdot = sinE*edot*(1.0 + satEph->e*cos(v)) / (sin(v)*(1.0-satEph->e*cosE));
	double udot = vdot + 2.0*(satEph->cus*cos2u - satEph->cuc*sin2u)*vdot;
	double rdot = satEph->A*satEph->e*sinE*edot + 2.0*(satEph->crs*cos2u - satEph->crc*sin2u)*vdot;
	double idotdot = satEph->idot + (satEph->cis*cos2u - satEph->cic*sin2u)*2*vdot;

	double vx_op = rdot*cos(u) - y_op*udot;
	double vy_op = rdot*sin(u) + x_op*udot;
	double omegadot = satEph->OMGd - CONST_OEDot;

	double tmpa = vx_op - y_op*cosi*omegadot;
	double tmpb = x_op*omegadot + vy_op*cosi - y_op*sini*idotdot;

	// Assign velocity states
	double state_xdot = tmpa * cos_omegak - tmpb * sin_omegak;
	satState->x_dot = state_xdot;
	double state_ydot = tmpa * sin_omegak + tmpb * cos_omegak;
	satState->y_dot = state_ydot;
	double state_zdot = vy_op*sini + y_op*cosi*idotdot;
	satState->z_dot = state_zdot;

	satState->delta_t_dot = clkd;


	// TODO: implement ionospheric corrections?

	return 0;
}


// satStatesArr: 	Pointer to an array of pointers to each batch of satStates
// satStatesCont: 	Pointer to a contiguous array of allocated satStates_t's to be broken up between the batches
// satStatePtrArr:	Pointer to a contiguous array of allocated state_t's
// initTime:		The GPS time that satStatesArr[0] will be started at
__global__ void
SPS_Init_Sat_Pos(dsp::utils::satState_t **satStatesArr, dsp::utils::satState_t *satStatesCont,
		dsp::utils::state_t<double> *satStatePtrArr, int initTime) {

    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int stride = blockDim.x * gridDim.x;

    unsigned int currBatch;
    unsigned int currSat;

    while (i < SATPOS_NUM_BATCHES) {
    	// Split the continuous array into an array of pointers
    	satStatesArr[i] = satStatesCont + i;

    	// Fill in the metadata
    	satStatesArr[i]->ready_d = 0;
    	satStatesArr[i]->replace_d = 1;
    	// Init at one full buffer set behind where you really want so Update_Sat_Pos forwards to the right location
    	satStatesArr[i]->startTime_d = initTime + (int)(((int)i-1-SATPOS_NUM_BATCHES) * SATPOS_BATCH_MIN * 60);

    	i += stride;
    }

    // Reset thread index for clean indexing of the next set
    i = blockIdx.x * blockDim.x + threadIdx.x;

    // Make sure the arrays have been assigned before continuing
    __syncthreads();

    while (i < SATPOS_NUM_BATCHES * CONST_PRN_MAX) {
        // Update index trackers
        currBatch = i / CONST_PRN_MAX;
        currSat = POSMOD(i, CONST_PRN_MAX);

        // Split the continuous array into an array of pointers
    	satStatesArr[currBatch]->satPosPtr_d[currSat] = satStatePtrArr + (i * SATPOS_BATCH_NUM_POS);
    	satStatesArr[currBatch]->valid_d[currSat] = 0;

    	i += stride;
    }

	return;
}





// Expecting a pointer to a single satState buffer -- the buffer to update
// Note that the time used to calculate each SatPos is the transmit time for that satellite
__global__ void
SPS_Update_Sat_Pos(dsp::utils::satState_t **satStates, dsp::utils::ephSet_t *navData, int navDataSize, int currBatch, dsp::utils::ephSet_t *debugEph) {

    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int stride = blockDim.x * gridDim.x;

    __shared__ int ephToUse[CONST_PRN_MAX];
    __shared__ int batchTime;

    int currPos;
    int currSat;

    int test = 0;	// Debug

    // Mark this buffer as in-progress for being updated
	if (i == 0) {
		satStates[currBatch]->ready_d = false;
		batchTime = satStates[currBatch]->startTime_d + (SATPOS_NUM_BATCHES * SATPOS_BATCH_MIN * 60);
	}

    __syncthreads();

	// First, get the eph data for these batches
	// (Batches should be small enough that it's reasonable for all elements to use the same eph)
	while (i < CONST_PRN_MAX) {

		// Initialize yourself first, since initializers are not allowed for shared memory
		ephToUse[i] = -1;

		// Check all eph's and find the one closest to the passed-in time if we have valid ephems for it
		// TODO: cleaner eph search that doesn't require the brute force search?
		//       or just stick with this to be thorough?
		for(int currIdx = 0; currIdx < navDataSize; currIdx++) {

			// Check if the ephems are valid for this PRN
			if (navData[currIdx].ephValid[i] &&
					fabs(navData[currIdx].ephToes - batchTime) <
					fabs(navData[ephToUse[i]].ephToes - batchTime)) {

				// Store the new best
				// Note: if there are no valid eph's for this PRN, this value will stay at the default of 0
				// The satellite positions for this PRN will not be updated, and the metadata will
				// make this PRN invalid
				ephToUse[i] = currIdx;
			}


			// Debug: force-choose the ephemerides
			if (navData[currIdx].ephToes == 417600 && navData[currIdx].ephValid[currSat]) {
				// Not elegant, but deep-copy this eph into the debugEph provided
				debugEph->ephToes 	= navData[currIdx].ephToes;
				debugEph->leaps 	= navData[currIdx].leaps;

				for (int i = 0; i < CONST_PRN_MAX; i++) {
					if (i < CONST_PRN_MAX) 	{ debugEph->eph[i] = navData[currIdx].eph[i]; debugEph->ephValid[i] = navData[currIdx].ephValid[i]; }
					if (i < 4) 				{ debugEph->utc_gps[i] = navData[currIdx].utc_gps[i]; }
					if (i < 8) 				{ debugEph->ion_gps[i] = navData[currIdx].ion_gps[i]; }
				}

				ephToUse[i] = currIdx;

				break;
			}
		}

		// After checking all the ephems, set the validity metadata
		satStates[currBatch]->valid_d[i] = navData[ephToUse[i]].ephValid[i];


		i += stride;
    }

	// Reset the threads for clean indexing of the next set
	i = blockIdx.x * blockDim.x + threadIdx.x;

    __syncthreads();

	// Update every element in the batch if needed
	while (i < SATPOS_BATCH_NUM_POS*CONST_PRN_MAX) {

		currSat = i / SATPOS_BATCH_NUM_POS;
		currPos = POSMOD(i, SATPOS_BATCH_NUM_POS);

		// Find the right eph data for this prn here
		if (satStates[currBatch]->valid_d[currSat]) {

			if (currBatch==0&&currSat==1&&currPos==0) {
				test = 1;
			}

			if(SPS_Get_Sat_Pos(satStates[currBatch]->satPosPtr_d[currSat] + currPos,
							   &(navData[ephToUse[currSat]].eph[currSat]),
							   batchTime + (SATPOS_BATCH_MIN*60)*(currPos / (double)(SATPOS_BATCH_NUM_POS)),
							   currBatch, currSat, currPos, navData[ephToUse[currSat]].ephToes ) == -1) {
				// If an error occurred processing this satPos, set this eph to invalid
				//cuPrintf("[SatPos] In SPS_Update_Sat_Pos, failed to converge on Kepler iteration at time %f.", ephTime);
				satStates[currBatch]->valid_d[currSat] = false;
			}
		}

		i += stride;
	}

	// Reset the thread index for clean indexing of the next set
	i = blockIdx.x * blockDim.x + threadIdx.x;


    __syncthreads();

    // Update metadata for this SatPos batch
    if (i == 0) {

    	satStates[currBatch]->ready_d = true;
    	satStates[currBatch]->replace_d = 0;
    	satStates[currBatch]->startTime_d = batchTime;
    }

	return;
}



dsp::SatPos::SatPos(){
	
    ModuleName = "SatPos";
    AllocateInputs(3);
    AllocateOutputs(4);

    Started = 0;

    /**
     * INPUT CONFIGURATION
     */
    // Configure inputs
    ConfigExpectedInput(0, "InitRXTime", 	DOUBLE_t, 		VALUE, 	VECTORLENGTH_ANY);
    ConfigExpectedInput(1, "InitEph",		UNDEFINED_t, 	EPHEMS, VECTORLENGTH_ANY);
    ConfigExpectedInput(2, "CurrRXTime",	DOUBLE_t,		VALUE,	VECTORLENGTH_ANY);
    
    InsertParam("OldEphTermTime", (void*)&oldEphTermTime, FLOAT_t, sizeof(double), sizeof(double));

    std::clog << "[" << ModuleName << "] Configured inputs" << std::endl;

    ConfigOutput(0, "SatPosPtrs", DOUBLE_t, VALUE, CUDA_DEVICE, VECTORLENGTH_ANY, NULL, 0); // Pointers to every buffer
    ConfigOutput(1, "SatPosArrSize", INT_t, VALUE, HOST, 1, NULL, 0); 						// Size of each buffer
    ConfigOutput(2, "SatPosBatchTime", 	DOUBLE_t, VALUE, HOST, 1, NULL, 0);					// The amount of time (s) covered by one batch
    ConfigOutput(3, "SatPosDebugEph", DOUBLE_t, VALUE, CUDA_DEVICE, 1, NULL, 0);			// Debug -- ephems that match PyGNSS
}


dsp::SatPos::~SatPos(){
	if (Started) Stop();
	delete [] ephSetArr;
    delete [] inputs;
    delete [] outputs;
    delete [] expectedInputs;
}


int 
dsp::SatPos::Start(void* cuFlowStream) {
	
	// Check module status and report accordingly
	if (Started) {
        std::clog << "[" << ModuleName << "] Start: Already Started." << std::endl;
        return 0;
    }
    std::clog << "[" << ModuleName << "] Starting ... " << std::flush;


	// Set the CUDA Stream for the GPU operations
    cuStream = (cudaStream_t*)cuFlowStream;


    // Get inputs
    initRXTime_a	= *((double*)(inputs[0]->Data));
    ephSetVec 		= *((std::vector<dsp::utils::ephSet_t>*)(inputs[1]->Data));
    ephSetSize	 	= ephSetVec.size();


    // Convert eph's from vector to array so cudaMemcpy can easily copy to device
    ephSetArr = new dsp::utils::ephSet_t[ephSetSize];
    for (int i = 0; i < ephSetSize; i++) { ephSetArr[i].copyInto(ephSetVec[i]); }


    // Initialize memory
    // This creates a giant block of contiguous memory for all satStates and all satellite positions
    // This needs to be broken down and have pointers assigned to appropriate locations -- to be done in SPS_Init_Sat_Pos
    size_t size = sizeof(dsp::utils::satState_t)*SATPOS_NUM_BATCHES;
	cuCheckMSt(cudaMalloc((void**)&satPosAll_d, size));

    size = sizeof(dsp::utils::state_t<double>)*SATPOS_NUM_BATCHES*SATPOS_BATCH_NUM_POS*CONST_PRN_MAX;
    cuCheckMSt(cudaMalloc((void**)&satPosPtrArrAll_d, size));

    size = sizeof(dsp::utils::satState_t*)*SATPOS_NUM_BATCHES;
    cuCheckMSt(cudaMalloc((void**)&satPosPtr_d, size));

    // Allocate space for every ephSet in navData
    size = sizeof(dsp::utils::ephSet_t)*ephSetSize;
    cuCheckMSt(cudaMalloc((void**)&ephSetPtr_d, size));
    cuCheckMSt(cudaMemcpyAsync(ephSetPtr_d, ephSetArr, size, cudaMemcpyHostToDevice));

    // Debug
    cuCheckMSt(cudaMalloc((void**)&debugEph_d, sizeof(dsp::utils::ephSet_t)));


    // Find the closest previous minute (in seconds)
    prevEphMinute = floor(initRXTime_a / 60) * 60;
    // And store the start time of the furthest buffer (in seconds)
    // -2 because one full buffer is initialized behind the buffer corresponding to prevEphMinute
    currFurthestBatchTime = prevEphMinute + (SATPOS_NUM_BATCHES-2) * SATPOS_BATCH_MIN * 60;
    // The 0th indexed batch is the oldest
    currOldestBatchIdx = 0;


    // Initialize the SatPos batches (compute these on the default stream since we need to block it)
    //SPS_Init_Sat_Pos<<<1, 128>>>(satPos_d, satPosAll_d, satPosPtrArrAll_d, prevEphMinute);
    SPS_Init_Sat_Pos<<<1, 128, 0, *cuStream>>>(satPosPtr_d, satPosAll_d, satPosPtrArrAll_d, prevEphMinute);
    for (int i = 0; i < SATPOS_NUM_BATCHES; i++) {
    	SPS_Update_Sat_Pos<<<1, 128, SATPOS_SHARED_SIZE, *cuStream>>>(satPosPtr_d, ephSetPtr_d, ephSetSize, i, debugEph_d);
    }


    // Now that space is allocated, assign the position buffer and size
    outputs[0].Data = satPosPtr_d;
    outputs[0].VectorLength = SATPOS_NUM_BATCHES;
    satPosArrSize = SATPOS_BATCH_NUM_POS;
    outputs[1].Data = &satPosArrSize;
    outputs[1].VectorLength = 1;
    satPosBatchTime = SATPOS_BATCH_MIN*60;
    outputs[2].Data = &satPosBatchTime;
    outputs[2].VectorLength = 1;

    // Debug outputs
    outputs[3].Data = debugEph_d;
    outputs[3].VectorLength = 1;

	
    // Make sure all GPU tasks have completed before continuing
    cuCheckMSt(cudaStreamSynchronize(*cuStream));
    cuCheckMSt(cudaDeviceSynchronize());


    // Create the batch calculation stream
    cuCheck(cudaStreamCreate(&satPosStream));

    // Signifies that the next call to update() will be the first after start()
    Started = 1;

    std::clog << "Started." << std::endl;
    return 0;
}


int 
dsp::SatPos::Stop() {
	
    int ret = 0;
    if (Started == 0) {
        std::clog << "[" << ModuleName << "] Stop: Wasn't Started." << std::endl;
        return 0;
    }
    std::clog << "[" << ModuleName << "] Stopping ... " << std::flush;

    // Delete streams
    cuCheckV(cudaStreamDestroy(satPosStream));

    // Free device memory
    cuCheckMSp(cudaFree((void*)satPosAll_d));
    cuCheckMSp(cudaFree((void*)satPosPtrArrAll_d));
    cuCheckMSp(cudaFree((void*)satPosPtr_d));
    cuCheckMSp(cudaFree((void*)ephSetPtr_d));
    cuCheckMSp(cudaFree((void*)debugEph_d));


    Started = 0;
    std::clog << "Stopped." << std::endl;

    return ret;
}
	

int
dsp::SatPos::Update(void* cuFlowStream) {

    if (Started == 0){
        std::cerr << "[" << ModuleName
                  << "] Error: Update() Failed due to SatPos not initialized"
                  << std::endl;
        return -1;
    }

    // If the rxTime is now on the last buffer, update the oldest buffer
    currRxTime = *((double*)(inputs[2]->Data));
    if (currRxTime > currFurthestBatchTime) {
        // Launch kernel to check if satellite positions need updated every iteration
        // (don't launch a kernel unless the stream is empty)
    	// Note: this scheme is a little goofy... semaphore-based locking to queue/unqueue update requests
    	// would be more robust
        if (cudaStreamQuery(satPosStream) == cudaSuccess) {
        	SPS_Update_Sat_Pos<<<1,128, SATPOS_SHARED_SIZE, satPosStream>>>(satPosPtr_d, ephSetPtr_d, ephSetSize, currOldestBatchIdx, debugEph_d);
        	std::clog << "[" << ModuleName << "] Updating SatPos Batch " << currOldestBatchIdx << std::endl;
        	// Move forward the track of the current oldest batch
        	currOldestBatchIdx = (currOldestBatchIdx + 1) % SATPOS_NUM_BATCHES;
        	currFurthestBatchTime += SATPOS_BATCH_MIN * 60;


        }
    }


	return 0;
}


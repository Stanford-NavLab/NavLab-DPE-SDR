
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
//#include "cuPrintf.cu"
#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <stdint.h>
#include "cuchanmgr.h"
#include "errorhandler.h"
#include "auxil.h"


__device__ __constant__ int 	CHM_LCA_d;
__device__ __constant__ double 	CHM_T_d;




/*
 *
 * DEVICE FUNCTIONS
 *
 */
__device__ double
CHM_Correct_Week_Crossover(double time) {
	if (time > 302400.0) 		{ return time-604800.0;	}
	else if (time < -302400.00)	{ return time+604800.0;	}
	else						{ return time;			}
}




// Light-weight, parallelized lat-lon calculator that returns values in radians
__device__ void
CHM_Dev_ECEF2LL_Rad(const double *posECEF, double *posLL) {

	// Compute lat
	double p = norm(2, posECEF);
	double theta = atan2(posECEF[2]*CONST_WGS84_A, p*CONST_WGS84_B);
	posLL[0] = atan2((posECEF[2] + pow(CONST_WGS84_EP,2) * CONST_WGS84_B * pow(sin(theta),3)),
			  (p - pow(CONST_WGS84_E,2) * CONST_WGS84_A * pow(cos(theta), 3)));

	// Compute lon
	posLL[1] = atan2(posECEF[1], posECEF[0]);

	return;
}


// Given latitude and longitude, return the elements of the ENU->ECEF rotation matrix
__device__ void
CHM_Dev_R_ENU2ECEF(const double *posLL, double *posENU) {

	double sinLat = sin(posLL[0]);
	double sinLon = sin(posLL[1]);
	double cosLat = cos(posLL[0]);
	double cosLon = cos(posLL[1]);

	posENU[0] = -sinLon;
	posENU[1] = -sinLat*cosLon;
	posENU[2] = cosLat*cosLon;
	posENU[3] = cosLon;
	posENU[4] = -sinLat*sinLon;
	posENU[5] = cosLat*sinLon;
	posENU[6] = 0.0;
	posENU[7] = cosLat;
	posENU[8] = sinLat;

	return;
}





// Function to be called within Update_Sat_Pos to prevent function saturation
// Note: this is hardcoded to compute GPS ephemerides only!
// See ephemeris.c in RTKLIB for potentially useful code for multi-constellation
//
// NOTE: IF THERE ARE ERRORS, "FMOD" is defined differently than np.mod! Check those lines!!!
//
__device__ int
CHM_Get_Sat_Pos(dsp::utils::state_t<double> *satState, dsp::utils::eph_t *satEph, double txTime) {

	// Corrected mean motion
	double n = sqrt(dsp::utils::MU_GPS/(CHM_CUBE(satEph->A))) + satEph->deln;

	// Compute satellite clock corrections
	double tc = CHM_Correct_Week_Crossover(txTime - satEph->tocs); 					// Without corrections
	double clkb = satEph->f2*tc*tc + satEph->f1*tc + satEph->f0 - satEph->tgd[0];	// Without relativistic correction
	double tk = CHM_Correct_Week_Crossover(txTime - clkb - satEph->toes);			// Without relativistic correction

	// Mean anomaly
	double E, M, f, dfdE, dE = 1;
	E = M = fmod((satEph->M0 + n*tk), CONST_2PI);
	// Eccentric anomaly
	for (int eccIdx = 0; eccIdx < CHM_MAX_ITER_KEPLER && abs(dE) > CHM_RTOL_KEPLER; eccIdx++) {
		f = M - E + satEph->e * sin(E);
		dfdE = -1.0 + satEph->e * cos(E);
		dE = -f / dfdE;
		E = fmod(E + dE, CONST_2PI);
		//if (abs(dE) < SATPOS_RTOL_KEPLER) { break; } // Break here if convergence is achieved!
	}
	if (abs(dE) > CHM_RTOL_KEPLER) { return -1;	}

	// Add in relativistic corrections and compute clock drift
	double dtr = CONST_F*(satEph->e)*(satEph->sqrt_A)*sin(E);
	tc = txTime - (clkb + dtr) - satEph->tocs;
	clkb = satEph->f2*tc*tc + satEph->f1*tc + satEph->f0 + dtr - satEph->tgd[0];
	double clkd = satEph->f1 + 2.0*satEph->f2*tc;


	// Recompute tk with relativisitic correction
	tk = CHM_Correct_Week_Crossover(txTime - clkb - satEph->toes);

	// Recompute mean anomaly with relativisitic correction
	dE = 1;
	E = M = fmod((satEph->M0 + n*tk), (CONST_2PI));
	// Eccentric anomaly
	for (int eccIdx = 0; eccIdx < CHM_MAX_ITER_KEPLER && abs(dE) > CHM_RTOL_KEPLER; eccIdx++) {
		f = M - E + satEph->e * sin(E);
		dfdE = -1.0 + satEph->e * cos(E);
		dE = -f / dfdE;
		E = fmod(E + dE, CONST_2PI);
		//if (abs(dE) < SATPOS_RTOL_KEPLER) { break; } // Break here if convergence is achieved!
	}
	if (abs(dE) > CHM_RTOL_KEPLER) { return -1;	}

	// Compute helpers
	double sinE = sin(E);
	double cosE = cos(E);


	// True anomaly
	double v = atan2(sqrt(1.0-satEph->e_sqr) * sinE / (1.0-satEph->e*cosE), (cosE-satEph->e)/(1.0-satEph->e*cosE));
	// Argument of latitude
	double u = fmod(v + satEph->omg, CONST_2PI); // Don't need mod here? (Computing for guarantee)

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






/*
 *
 * DEVICE KERNELS
 *
 */

/** \brief Compute the satellite state given the current estimate of the received signal
 *
 * The initial update -- propagates forward from scalar handoff
 *
 * \param navData		Ptr to all ephemerides we have
 * \param navDataSize	The number of ephemerides in navData
 * \param ephToUse		Output: the ephemerides chosen (closest in time to the current state estimate)
 * \param satStates		Output: The states of the tracked satellites according to the channel parameters
 * \param codePhase		Ptr to the current estimate of the code phase of each channel
 * \param cpElapsed		The number of code periods for this sample set since the start of tracking
 * \param cpRef			The number of code periods elapsed since the start of tracking for the reference TOW
 * \param cpRefTOW		The TOW at the reference code period
 * \param PRNs			The satellite channels being tracked
 * \param numChan		The number of channels being tracked
 * \param txTime		Output: When each satellite sent the transmission that is now the sample set (computed by the channel parameters)
 *
 */
__global__ void
CHM_ComputeSatStates(dsp::utils::ephSet_t *navData, int navDataSize, int *ephToUse,
		dsp::utils::state_t<double> *satStates,
		double *codePhase,
		int *cpElapsed, int *cpRef, int *cpRefTOW,
		uint8_t *PRNs, int numChan,
		double *txTime) {

	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int stride = blockDim.x * gridDim.x;

	int test = 0;


    // only compute if within the range to estimate
    while (i < numChan) {

	    // COMPUTE TXTIME
		txTime[i] = cpRefTOW[i] +
					((cpElapsed[i] - cpRef[i]) * CONST_T_CA) +
					(codePhase[i] / CONST_F_CA);



	    // SATELLITE COMPUTATION

    	// Save the PRN number for this index
    	int currPrn = PRNs[i]-1;

		// Reset eph to the unchosen index
		ephToUse[i] = -1;


		// Check all eph's and find the one closest to the passed-in time if we have valid ephems for it
		// TODO: cleaner eph search that doesn't require the brute force search?
		//       or just stick with this to be thorough?
		for(int currIdx = 0; currIdx < navDataSize; currIdx++) {

			// Check if the ephems are valid for this PRN
			if (navData[currIdx].ephValid[currPrn]) {

				// If no eph chosen and this eph is valid, select this eph
				if (ephToUse[i] == -1) {
					ephToUse[i] = currIdx;
				}
				// Update the best eph if it's closer in time to the transmit time
				else if (fabs(navData[currIdx].ephToes  - txTime[i]) <
						 fabs(navData[ephToUse[i]].ephToes - txTime[i])) {
					ephToUse[i] = currIdx;
				}
			}

		}

		// Make sure that an ephemeris actually got selected
		if (ephToUse[i] == -1) {
			// return error
			// FIXME: Remove channel
			test = 1;
		}

		CHM_Get_Sat_Pos(&(satStates[i]), &(navData[ephToUse[i]].eph[currPrn]), txTime[i]);


    	i += stride;
    }
}




/** \brief The standard update: Performs MeasUpdate -> TimeUpdate -> txTime -> satState
 *
 * Propagates from EKF measurement; to be run in update
 *
 * \param navData			Ptr to all ephemerides we have
 * \param navDataSize		The number of ephemerides in navData
 * \param ephToUse			The ephemerides chosen (closest in time to the current state estimate)
 * \param centerPt			The current best estimate on the receiver's state
 * \param satStates			Output: the states of the tracked satellites according to the channel parameters
 * \param codePhaseStart	Ptr to the estimate of the code phase of each channel at the beginning of the sample set
 * \param codePhaseEnd		Ptr to the estimate of the code phase of each channel at the end of the sample set
 * \param codeFreq			Ptr to the estimate of code frequency (under Doppler effects)
 * \param carrPhaseStart	Ptr to the estimate of the carrier phase of each channel at the beginning of the sample set
 * \param carrPhaseEnd		Ptr to the estimate of the carrier phase of each channel at the end of the sample set
 * \param carrFreq			Ptr to the estimate of the carrier frequency (under Doppler effects)
 * \param dopplerSign		The sign convention for Doppler frequency
 * \param cpElapsedStart	The number of code periods for the start of this sample set since the start of tracking
 * \param cpElapsedEnd		The number of code periods for the end of this sample set since the start of tracking
 * \param cpRef				The number of code periods elapsed since the start of tracking for the reference TOW
 * \param cpRefTOW			The TOW at the reference code period
 * \param PRNs				The satellite channels being tracked
 * \param numChan			The number of channels being tracked
 * \param rxTime			The current time of the receiver's internal clock (no delta-t applied)
 * \param txTime			Output: When each satellite sent the transmission that is now the sample set (computed by the channel parameters)
 * \param T					The length in seconds of one sample set
 *
 */
__global__ void
CHM_PropagateChannels(dsp::utils::ephSet_t *navData, int navDataSize, int *ephToUse,
		const double *centerPt, dsp::utils::state_t<double> *satStates,
		double *codePhaseStart, double *codePhaseEnd, double *codeFreq,
		double *carrPhaseStart, double *carrPhaseEnd, double *carrFreq, int *dopplerSign,
		int *cpElapsedStart, int *cpElapsedEnd, int *cpRef, int *cpRefTOW,
		uint8_t *PRNs, int numChan,
		double rxTime, double *txTime, double T)
{
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int stride = blockDim.x * gridDim.x;

	int test = 0;

    // Back-calculation values for finding scores
    int currPrn;
    double currPtVelECI[4];
    double satPosTransmitTime;
    dsp::utils::state_t<double> currPtSatState;
    double bc_los[3];
    double bc_range;
    double bc_losrangerate;
    double bc_fi;
    double bc_rc0;
    double bc_fc;
	double bc_pseudorange;
	double bc_txTime;
	double bc_codeFracDiff;
	double bc_rc;
    double cos_tau_OEDot;
    double sin_tau_OEDot;



    // Since kernels are launched in multiples of 32 threads for efficiency,
    // only compute if within the range to estimate
    while (i < numChan) {


    	// MEASUREMENT UPDATE

		// For this satellite position, compute the transmit time (TOF) to the candidate point
		satPosTransmitTime = rxTime - (txTime[i] + (centerPt[3]/(double)CONST_C)) + satStates[i].delta_t;

		// Convert the satellite and coordinate positions to ECI to add
		cos_tau_OEDot = cos(-CONST_OEDot*satPosTransmitTime);
		sin_tau_OEDot = sin(-CONST_OEDot*satPosTransmitTime);


		// Rotate the satellite position over the earth for this point's pvt state
		currPtSatState.x 			= cos_tau_OEDot * satStates[i].x
									- sin_tau_OEDot * satStates[i].y;
		currPtSatState.y 			= sin_tau_OEDot * satStates[i].x
									+ cos_tau_OEDot * satStates[i].y;
		currPtSatState.z 			= satStates[i].z;
		currPtSatState.delta_t 		= satStates[i].delta_t;

		currPtSatState.x_dot		= cos_tau_OEDot * satStates[i].x_dot
									- sin_tau_OEDot * satStates[i].y_dot
									- CONST_OEDot * sin_tau_OEDot * satStates[i].x
									- CONST_OEDot * cos_tau_OEDot * satStates[i].y;
		currPtSatState.y_dot		= sin_tau_OEDot * satStates[i].x_dot
									+ cos_tau_OEDot * satStates[i].y_dot
									+ CONST_OEDot * cos_tau_OEDot * satStates[i].x
									- CONST_OEDot * sin_tau_OEDot * satStates[i].y;
		currPtSatState.z_dot		= satStates[i].z_dot;
		currPtSatState.delta_t_dot	= satStates[i].delta_t_dot;


		// Also need to rotate the point's velocity into the inertial frame
		currPtVelECI[0] = centerPt[4] - CONST_OEDot*centerPt[1];
		currPtVelECI[1] = centerPt[5] + CONST_OEDot*centerPt[0];
		currPtVelECI[2] = centerPt[6];
		currPtVelECI[3] = centerPt[7];



		// Back-calculate the channel parameters

		// Find carrier frequency
		bc_los[0] = currPtSatState.x - centerPt[0];
		bc_los[1] = currPtSatState.y - centerPt[1];
		bc_los[2] = currPtSatState.z - centerPt[2];
		bc_range = norm(3, bc_los); // Might as well find range this way since the unit vector is needed.
		bc_losrangerate = 	((bc_los[0]/bc_range)*(currPtVelECI[0]-currPtSatState.x_dot)) +
							((bc_los[1]/bc_range)*(currPtVelECI[1]-currPtSatState.y_dot)) +
							((bc_los[2]/bc_range)*(currPtVelECI[2]-currPtSatState.z_dot));
		bc_fi = CONST_F_L1 * ((bc_losrangerate - currPtVelECI[3])/CONST_C + currPtSatState.delta_t_dot) / (*dopplerSign);


		// Find the code frequency elapsed since last timestep
		bc_pseudorange = bc_range - CONST_C * currPtSatState.delta_t + centerPt[3];
		bc_txTime = rxTime - bc_pseudorange/CONST_C;
		bc_codeFracDiff = bc_txTime - cpRefTOW[i] - ((cpElapsedEnd[i] - cpRef[i]) * CONST_T_CA);
		bc_rc = bc_codeFracDiff * CONST_F_CA;
		// The back-calculated code frequency is the updated frequency for the channel
		bc_fc = CONST_F_CA + (*dopplerSign * CONST_F_CA/CONST_F_L1) * bc_fi + (bc_rc - codePhaseEnd[i]) / CHM_T_d;



		// Something's very wrong if the codeFreq is shifted by more than 10kHz from F_CA
		if (abs(bc_fc-CONST_F_CA) > 10000) {
			// TODO: Throw an error if we get here
			test = 3;
		}


		// Update channels using the calculated results
		carrFreq[i] = bc_fi;
		codeFreq[i] = bc_fc;



    	// TIME UPDATE

    	// ***ENHANCED*** time update -- uses the back-calculated code phase instead
    	// Progress the channel params by one sample set (not ephem dependent; can still do regardless)
    	double temp1 = floor((codeFreq[i]*CHM_T_d + codePhaseEnd[i])/CHM_LCA_d);
    	if (temp1 > 30 || temp1 < 10) {
    		test = -1;
    	}
    	double cpElapsedPred = cpElapsedEnd[i] + temp1;
    	double temp2 = fmod((double)(codeFreq[i]*CHM_T_d + codePhaseEnd[i]), (double)CHM_LCA_d);
    	if (temp2 < 0.0) {
    		temp2 += (double)CHM_LCA_d;
    	}
    	double codePhasePred = temp2;

	    // COMPUTE TXTIME
    	// Reset the PRN flag
		double txTimePred = cpRefTOW[i] +
					((cpElapsedPred - cpRef[i]) * CONST_T_CA) +
					(codePhasePred / CONST_F_CA);



	    // SATELLITE COMPUTATION

    	// Save the PRN number for this index
    	int currPrn = PRNs[i]-1;

		// Reset eph to the unchosen index
		ephToUse[i] = -1;


		// Check all eph's and find the one closest to the passed-in time if we have valid ephems for it
		// TODO: cleaner eph search that doesn't require the brute force search?
		//       or just stick with this to be thorough?
		for(int currIdx = 0; currIdx < navDataSize; currIdx++) {

			// Check if the ephems are valid for this PRN
			if (navData[currIdx].ephValid[currPrn]) {

				// If no eph chosen and this eph is valid, select this eph
				if (ephToUse[i] == -1) {
					ephToUse[i] = currIdx;
				}
				// Update the best eph if it's closer in time to the transmit time
				else if (fabs(navData[currIdx].ephToes  - txTimePred) <
						 fabs(navData[ephToUse[i]].ephToes - txTimePred)) {
					ephToUse[i] = currIdx;
				}
			}
		}

		// Make sure that an ephemeris actually got selected
		if (ephToUse[i] == -1) {
			// return error
			// FIXME: Remove channel
			test = 1;
		}

		dsp::utils::state_t<double> satStatePred;
		CHM_Get_Sat_Pos(&satStatePred, &(navData[ephToUse[i]].eph[currPrn]), txTimePred);


		// (Internally) advance the rxTime by T and recompute the code phase
		// For this satellite position, compute the transmit time (TOF) to the candidate point
		satPosTransmitTime = rxTime + T - (txTimePred + (centerPt[3]/(double)CONST_C)) + satStatePred.delta_t;

		// Convert the satellite and coordinate positions to ECI to add
		cos_tau_OEDot = cos(-CONST_OEDot*satPosTransmitTime);
		sin_tau_OEDot = sin(-CONST_OEDot*satPosTransmitTime);

		// Rotate the satellite position over the earth for this point's pvt state
		currPtSatState.x 			= cos_tau_OEDot * satStatePred.x
									- sin_tau_OEDot * satStatePred.y;
		currPtSatState.y 			= sin_tau_OEDot * satStatePred.x
									+ cos_tau_OEDot * satStatePred.y;
		currPtSatState.z 			= satStatePred.z;
		currPtSatState.delta_t 		= satStatePred.delta_t;

		// Back-calculate the channel parameters
		// Find carrier frequency
		bc_los[0] = currPtSatState.x - centerPt[0];
		bc_los[1] = currPtSatState.y - centerPt[1];
		bc_los[2] = currPtSatState.z - centerPt[2];
		bc_range = norm(3, bc_los); // Might as well find range this way since the unit vector is needed.

		bc_pseudorange = bc_range - CONST_C * currPtSatState.delta_t + centerPt[3];
		bc_txTime = rxTime + T - bc_pseudorange/CONST_C;
		bc_codeFracDiff = bc_txTime - cpRefTOW[i] - ((cpElapsedEnd[i] - cpRef[i]) * CONST_T_CA);
		bc_rc = bc_codeFracDiff * CONST_F_CA;

		double bc_rc0_check = bc_rc - codePhaseEnd[i];

		// Find the code frequency elapsed since last timestep
		bc_rc0 = CONST_F_CA * (satPosTransmitTime - (bc_range / CONST_C));



    	// Debug comparison
    	double temp1old = floor((codeFreq[i]*CHM_T_d + codePhaseEnd[i])/CHM_LCA_d);
    	if (temp1old > 30 || temp1old < 10) {
    		test = -1;
    	}
    	double temp2old = fmod((double)(codeFreq[i]*CHM_T_d + codePhaseEnd[i]), (double)CHM_LCA_d);
    	if (temp2old < 0.0) {
    		temp2old += (double)CHM_LCA_d;
    	}



    	// Progress the channel params by one sample set (not ephem dependent; can still do regardless)
		cpElapsedStart[i] = cpElapsedEnd[i];
    	temp1 = floor(bc_rc/CHM_LCA_d);
		if (temp1 > 30 || temp1 < 10) {
			test = -1;
		}

		codePhaseStart[i] = codePhaseEnd[i];
		temp2 = fmod(bc_rc, (double)CHM_LCA_d);
		if (temp2 < 0.0) {
			temp2 += (double)CHM_LCA_d;
		}
		cpElapsedEnd[i] += temp1;
		codePhaseEnd[i] = temp2;

		carrPhaseStart[i] = carrPhaseEnd[i];
		double temp3 = fmod((double)(carrFreq[i]*CHM_T_d + carrPhaseEnd[i]), (double)1.0);
		if (temp3 < 0.0) {
			temp3 += 1.0;
		}
		carrPhaseEnd[i] = temp3;

		// End ENHANCED




	    // COMPUTE TXTIME
    	// Reset the PRN flag
		txTime[i] = cpRefTOW[i] +
					((cpElapsedEnd[i] - cpRef[i]) * CONST_T_CA) +
					(codePhaseEnd[i] / CONST_F_CA);




	    // SATELLITE COMPUTATION

    	// Save the PRN number for this index
    	currPrn = PRNs[i]-1;

		CHM_Get_Sat_Pos(&(satStates[i]), &(navData[ephToUse[i]].eph[currPrn]), txTime[i]);

    	i += stride;
    }

	return;
}





/** \brief The initial update: Performs TimeUpdate -> txTime -> satState
 *
 * Propagates from EKF measurement; to be run in update
 *
 * \param navData			Ptr to all ephemerides we have
 * \param navDataSize		The number of ephemerides in navData
 * \param ephToUse			The ephemerides chosen (closest in time to the current state estimate)
 * \param centerPt			The current best estimate on the receiver's state
 * \param satStates			Output: the states of the tracked satellites according to the channel parameters
 * \param codePhaseStart	Ptr to the estimate of the code phase of each channel at the beginning of the sample set
 * \param codePhaseEnd		Ptr to the estimate of the code phase of each channel at the end of the sample set
 * \param codeFreq			Ptr to the estimate of code frequency (under Doppler effects)
 * \param carrPhaseStart	Ptr to the estimate of the carrier phase of each channel at the beginning of the sample set
 * \param carrPhaseEnd		Ptr to the estimate of the carrier phase of each channel at the end of the sample set
 * \param carrFreq			Ptr to the estimate of the carrier frequency (under Doppler effects)
 * \param dopplerSign		The sign convention for Doppler frequency
 * \param cpElapsedStart	The number of code periods for the start of this sample set since the start of tracking
 * \param cpElapsedEnd		The number of code periods for the end of this sample set since the start of tracking
 * \param cpRef				The number of code periods elapsed since the start of tracking for the reference TOW
 * \param cpRefTOW			The TOW at the reference code period
 * \param PRNs				The satellite channels being tracked
 * \param numChan			The number of channels being tracked
 * \param rxTime			The current time of the receiver's internal clock (no delta-t applied)
 * \param txTime			Output: When each satellite sent the transmission that is now the sample set (computed by the channel parameters)
 * \param T					The length in seconds of one sample set
 *
 */
__global__ void
CHM_TimeUpdateChannels(dsp::utils::ephSet_t *navData, int navDataSize, int *ephToUse,
		const double *centerPt, dsp::utils::state_t<double> *satStates,
		double *codePhaseStart, double *codePhaseEnd, double *codeFreq,
		double *carrPhaseStart, double *carrPhaseEnd, double *carrFreq, int *dopplerSign,
		int *cpElapsedStart, int *cpElapsedEnd, int *cpRef, int *cpRefTOW,
		uint8_t *PRNs, int numChan,
		double rxTime, double *txTime, double T)
{
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int stride = blockDim.x * gridDim.x;

	int test = 0;

    // Back-calculation values for finding scores
    int currPrn;
    double currPtVelECI[4];
    double satPosTransmitTime;
    dsp::utils::state_t<double> currPtSatState;
    double bc_los[3];
    double bc_range;
    double bc_losrangerate;
    double bc_fi;
    double bc_rc0;
    double bc_fc;
    double cos_tau_OEDot;
    double sin_tau_OEDot;


    // Since kernels are launched in multiples of 32 threads for efficiency,
    // only compute if within the range to estimate
    while (i < numChan) {


    	// TIME UPDATE
    	// ***ENHANCED*** time update -- uses the back-calculated code phase instead
    	// Progress the channel params by one sample set (not ephem dependent; can still do regardless)
    	//cpElapsed[i]   += floor((BCS_S_d * (codeFreq[i]/BCS_FS_d) + codePhase[i])/BCS_LCA_d);
    	double temp1 = floor((codeFreq[i]*CHM_T_d + codePhaseEnd[i])/CHM_LCA_d);
    	if (temp1 > 30 || temp1 < 10) {
    		test = -1;
    	}
    	double cpElapsedPred = cpElapsedEnd[i] + temp1;
    	//codePhase[i] 	= fmod((double)(codePhase[i] + codeFreq[i]*BCS_T_d), (double)BCS_LCA_d);
    	double temp2 = fmod((double)(codeFreq[i]*CHM_T_d + codePhaseEnd[i]), (double)CHM_LCA_d);
    	if (temp2 < 0.0) {
    		temp2 += (double)CHM_LCA_d;
    	}
    	double codePhasePred = temp2;



	    // COMPUTE TXTIME
    	// Reset the PRN flag
		double txTimePred = cpRefTOW[i] +
					((cpElapsedPred - cpRef[i]) * CONST_T_CA) +
					(codePhasePred / CONST_F_CA);



	    // SATELLITE COMPUTATION

    	// Save the PRN number for this index
    	int currPrn = PRNs[i]-1;

		// Reset eph to the unchosen index
		ephToUse[i] = -1;


		// Check all eph's and find the one closest to the passed-in time if we have valid ephems for it
		// TODO: cleaner eph search that doesn't require the brute force search?
		//       or just stick with this to be thorough?
		for(int currIdx = 0; currIdx < navDataSize; currIdx++) {

			// Check if the ephems are valid for this PRN
			if (navData[currIdx].ephValid[currPrn]) {

				// If no eph chosen and this eph is valid, select this eph
				if (ephToUse[i] == -1) {
					ephToUse[i] = currIdx;
				}
				// Update the best eph if it's closer in time to the transmit time
				else if (fabs(navData[currIdx].ephToes  - txTimePred) <
						 fabs(navData[ephToUse[i]].ephToes - txTimePred)) {
					ephToUse[i] = currIdx;
				}
			}
		}

		// Make sure that an ephemeris actually got selected
		if (ephToUse[i] == -1) {
			// TODO: Return error
			test = 1;
		}

		dsp::utils::state_t<double> satStatePred;
		CHM_Get_Sat_Pos(&satStatePred, &(navData[ephToUse[i]].eph[currPrn]), txTimePred);


		// (Internally) advance the rxTime by T and recompute the code phase
		// For this satellite position, compute the transmit time (TOF) to the candidate point
		satPosTransmitTime = rxTime + T - (txTimePred + (centerPt[3]/(double)CONST_C)) + satStatePred.delta_t;

		// Convert the satellite and coordinate positions to ECI to add
		cos_tau_OEDot = cos(-CONST_OEDot*satPosTransmitTime);
		sin_tau_OEDot = sin(-CONST_OEDot*satPosTransmitTime);

		// Rotate the satellite position over the earth for this point's pvt state
		currPtSatState.x 			= cos_tau_OEDot * satStatePred.x
									- sin_tau_OEDot * satStatePred.y;
		currPtSatState.y 			= sin_tau_OEDot * satStatePred.x
									+ cos_tau_OEDot * satStatePred.y;
		currPtSatState.z 			= satStatePred.z;
		currPtSatState.delta_t 		= satStatePred.delta_t;

		// Back-calculate the channel parameters
		// Find carrier frequency
		bc_los[0] = currPtSatState.x - centerPt[0];
		bc_los[1] = currPtSatState.y - centerPt[1];
		bc_los[2] = currPtSatState.z - centerPt[2];
		bc_range = norm(3, bc_los); // Might as well find range this way since the unit vector is needed.

		double bc_pseudorange = bc_range - CONST_C * currPtSatState.delta_t + centerPt[3];

		double bc_txTime = rxTime + T - bc_pseudorange/CONST_C;

		// Not using the predicted values here because we want to see how much the code phase moved since last timestep
		// based on the predicted satellite location and newly-predicted receiver location.
		// This would be most accurate if done iteratively, as the new code phase estimate would update the satellite position.
		// But, after re-applying the code phase estimate to the satellite position and calculating again, the bc-pseudorange
		// shouldn't change much, because the satellites move at 4m/ms and we expect code chips of difference with bc_rc (1ms/1023)
		double bc_codeFracDiff = bc_txTime - cpRefTOW[i] - ((cpElapsedEnd[i] - cpRef[i]) * CONST_T_CA);

		double bc_rc = bc_codeFracDiff * CONST_F_CA;

		double bc_rc0_check = bc_rc - codePhaseEnd[i];
		double modCheck = fmod(bc_rc0_check, (double)CHM_LCA_d);

		// Find the code frequency elapsed since last timestep
		bc_rc0 = CONST_F_CA * (satPosTransmitTime - (bc_range / CONST_C));

		// Advance the code phase
    	double temp2old = fmod((double)(codeFreq[i]*CHM_T_d + codePhaseEnd[i]), (double)CHM_LCA_d);
    	if (temp2old < 0.0) {
    		temp2old += (double)CHM_LCA_d;
    	}

    	// Progress the channel params by one sample set (not ephem dependent; can still do regardless)
		cpElapsedStart[i] = cpElapsedEnd[i];
    	temp1 = floor(bc_rc/CHM_LCA_d);
		if (temp1 > 30 || temp1 < 10) {
			test = -1;
		}
		codePhaseStart[i] = codePhaseEnd[i];
		temp2 = fmod(bc_rc, (double)CHM_LCA_d);
		if (temp2 < 0.0) {
			temp2 += (double)CHM_LCA_d;
		}
		cpElapsedEnd[i] += temp1;
		codePhaseEnd[i] = temp2;

		carrPhaseStart[i] = carrPhaseEnd[i];
		double temp3 = fmod((double)(carrFreq[i]*CHM_T_d + carrPhaseEnd[i]), (double)1.0);
		if (temp3 < 0.0) {
			temp3 += 1.0;
		}
		carrPhaseEnd[i] = temp3;


	    // COMPUTE TXTIME
    	// Reset the PRN flag
		txTime[i] = cpRefTOW[i] +
					((cpElapsedEnd[i] - cpRef[i]) * CONST_T_CA) +
					(codePhaseEnd[i] / CONST_F_CA);



	    // SATELLITE COMPUTATION

    	// Save the PRN number for this index
    	currPrn = PRNs[i]-1;

		CHM_Get_Sat_Pos(&(satStates[i]), &(navData[ephToUse[i]].eph[currPrn]), txTime[i]);

    	i += stride;
    }

	return;
}










/** \brief Computes values needed by BCM grid evaluation -- ENU->ECEF rotation matrix and the satellite states
 *
 * \param rxTime			The current time of the receiver's internal clock (no delta-t applied)
 * \param txTime			When each satellite sent the transmission that is now the sample set (computed by the channel parameters)
 * \param centerPt			The current best estimate on the receiver's state
 * \param satStates			The states of the tracked satellites according to the channel parameters
 * \param numChan			The number of channels being tracked
 * \param timeGrid			The specific offsets from centerPt being evaluated on the grid
 * \param timeGridDim		The number of time offsets being evaluated
 * \param batchSatStates	Output: The specific states of the satellites at all the offsets on the grid
 * \param enu2ecefMat		Output: The rotation matrix to convert the centerPt-ENU frame to ECEF
 *
 */
__global__ void
CHM_GridPrep(double rxTime, double *txTime, double *centerPt, dsp::utils::state_t<double> *satStates,
		int numChan, double *timeGrid, int timeGridDim,
		dsp::utils::state_t<double> *batchSatStates, double *enu2ecefMat) {

	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int stride = blockDim.x * gridDim.x;

	int currPt; 	// The point on the grid being examined
	int currChan; 	// The channel this thread is processing

	double satPosTransmitTime;
	double cos_tau_OEDot;
	double sin_tau_OEDot;


	while(i < numChan * timeGridDim + 1) {

		// Update index tracker
		currPt = POSMOD(i, timeGridDim);
		currChan = i / timeGridDim;


		// If the last thread, compute the ENU->ECEF rotation matrix
		// This looks like it goes ECEF->ENU, but ENU->ECEF is the transpose of ECEF->ENU,
		// and ENU->ECEF gives the transformation from grid-space to satellite frame space,
		// so that's the rotation we actually want.
		if (i == numChan * timeGridDim) {
			double posLL[2];
			CHM_Dev_ECEF2LL_Rad(centerPt, posLL);
			CHM_Dev_R_ENU2ECEF(posLL, enu2ecefMat);
		}
		// Otherwise, find the satellite state you're responsible for rotating
		// We know the satellite's location from the channel parameters -- it's a function of the received signal.
		// However, the state we compute from the channel parameters needs to consider the rotation of the Earth
		// for EACH time being evaluated on the grid. So, the multiple satellite states being computed here are the same
		// location in different ECEF frame-times.
		else {
    		// For this satellite position, compute the transmit time (TOF) to the candidate point
    		satPosTransmitTime = rxTime - (txTime[currChan] + ((timeGrid[currPt] + centerPt[3])/(double)CONST_C)) + satStates[currChan].delta_t; // More conceptually accurate

    		// Convert the satellite and coordinate positions to ECI to add
			cos_tau_OEDot = cos(-CONST_OEDot*satPosTransmitTime);
			sin_tau_OEDot = sin(-CONST_OEDot*satPosTransmitTime);

			// batchSatStates contains the states for all grid points for the first channel, then for the second channel, ...
			// Rotate the satellite position over the earth for this point's pvt state
			batchSatStates[i].x 			= cos_tau_OEDot * satStates[currChan].x
											- sin_tau_OEDot * satStates[currChan].y;
			batchSatStates[i].y 			= sin_tau_OEDot * satStates[currChan].x
											+ cos_tau_OEDot * satStates[currChan].y;
			batchSatStates[i].z 			= satStates[currChan].z;
			batchSatStates[i].delta_t 		= satStates[currChan].delta_t;

			batchSatStates[i].x_dot			= cos_tau_OEDot * satStates[currChan].x_dot
    										- sin_tau_OEDot * satStates[currChan].y_dot
    										- CONST_OEDot * sin_tau_OEDot * satStates[currChan].x
    										- CONST_OEDot * cos_tau_OEDot * satStates[currChan].y;
			batchSatStates[i].y_dot			= sin_tau_OEDot * satStates[currChan].x_dot
    										+ cos_tau_OEDot * satStates[currChan].y_dot
    										+ CONST_OEDot * cos_tau_OEDot * satStates[currChan].x
    										- CONST_OEDot * sin_tau_OEDot * satStates[currChan].y;
			batchSatStates[i].z_dot			= satStates[currChan].z_dot;
			batchSatStates[i].delta_t_dot 	= satStates[currChan].delta_t_dot;

		}


		i += stride;
	}
}






dsp::cuChanMgr::cuChanMgr(){

    ModuleName = "cuChanMgr";
    AllocateInputs(15);
    AllocateOutputs(18);

    Started = 0;

    /**
     * INPUT CONFIGURATION
     */
    // Configure inputs
    ConfigExpectedInput(0, "InitEph",					UNDEFINED_t, 	EPHEMS, 		VECTORLENGTH_ANY);
    ConfigExpectedInput(1, "InitPRN", 					CHAR_t, 		VALUE, 			VECTORLENGTH_ANY);
    ConfigExpectedInput(2, "InitCodePhase", 			DOUBLE_t, 		VALUE, 			VECTORLENGTH_ANY);
    ConfigExpectedInput(3, "InitCarrierPhase", 			DOUBLE_t, 		VALUE, 			VECTORLENGTH_ANY);
    ConfigExpectedInput(4, "InitCodeFrequency", 		DOUBLE_t, 		FREQUENCY_HZ, 	VECTORLENGTH_ANY);
    ConfigExpectedInput(5, "InitCarrierFrequency", 		DOUBLE_t, 		FREQUENCY_HZ, 	VECTORLENGTH_ANY);
    ConfigExpectedInput(6, "InitElapsedCodePeriods", 	INT_t, 			VALUE, 			VECTORLENGTH_ANY);
    ConfigExpectedInput(7, "InitReferenceCodePeriods", 	INT_t, 			VALUE, 			VECTORLENGTH_ANY);
    ConfigExpectedInput(8, "InitCPRefTOW", 				INT_t, 			VALUE, 			VECTORLENGTH_ANY);
    ConfigExpectedInput(9, "InitRXTime", 				DOUBLE_t, 		VALUE, 			VECTORLENGTH_ANY);
    ConfigExpectedInput(10, "xCurrk1k1",				DOUBLE_t,		STATE,			VECTORLENGTH_ANY);
    ConfigExpectedInput(11, "SampleLength", 			DOUBLE_t, 		VALUE, 			1);
    ConfigExpectedInput(12, "xCurrkk1", 				DOUBLE_t,		STATE,			VECTORLENGTH_ANY);
    ConfigExpectedInput(13, "TimeGrid",					DOUBLE_t,		VALUE,			VECTORLENGTH_ANY);



    InsertParam("DopplerSign", 		(void*)&dopplerSign, 			INT_t, 		sizeof(int), 	sizeof(int));


    std::clog << "[" << ModuleName << "] Configured inputs" << std::endl;


    // Note: Phase and cpEla params are saved wrt both the start AND end of the sample set!
    // 		 BatchCorrScores wants references wrt the start so indices can be counted forwards
    //       State estimation is done wrt the end, since we have that full set of samples
    //       Frequency measurements are wrt the end, though this really just means they include the most recent measurement
    //       This is consistent with PyGNSS:
    //       rxTime+=T -> TimeUpdateState -> BatchCorr -> TimeUpdateChannels -> ManifoldEstimation -> MeasUpdateState -> MeasUpdateChannels
    //       PyGNSS does it that way so BatchCorr uses start-referenced params and ManifoldEstimation uses end-referenced params.
    //		 Though, doing the param updating all at once and denoting which is start-referenced and end-referenced is clearer.
    ConfigOutput(0, "rxTime", 				DOUBLE_t,	VALUE,			HOST,		 1,					NULL, 0);
    ConfigOutput(1, "txTime",     			DOUBLE_t,   VALUE,      	CUDA_DEVICE, VECTORLENGTH_ANY, 	NULL, 0);
    ConfigOutput(2, "CodePhaseStart",		DOUBLE_t,	VALUE,			CUDA_DEVICE, VECTORLENGTH_ANY,	NULL, 0);
    ConfigOutput(3, "CarrierPhaseStart",	DOUBLE_t,	VALUE,			CUDA_DEVICE, VECTORLENGTH_ANY,	NULL, 0);
    ConfigOutput(4, "CodePhaseEnd",			DOUBLE_t,	VALUE,			CUDA_DEVICE, VECTORLENGTH_ANY,	NULL, 0);
    ConfigOutput(5, "CarrierPhaseEnd",		DOUBLE_t,	VALUE,			CUDA_DEVICE, VECTORLENGTH_ANY,	NULL, 0);
    ConfigOutput(6, "CodeFrequency",		DOUBLE_t,	FREQUENCY_HZ,	CUDA_DEVICE, VECTORLENGTH_ANY,	NULL, 0);
    ConfigOutput(7, "CarrierFrequency",		DOUBLE_t,	FREQUENCY_HZ,	CUDA_DEVICE, VECTORLENGTH_ANY,	NULL, 0);
    ConfigOutput(8, "SatStates", 			DOUBLE_t,	STATE,			CUDA_DEVICE, VECTORLENGTH_ANY,	NULL, 0);
    ConfigOutput(9, "DopplerSign",			INT_t, 		VALUE,			CUDA_DEVICE, 1, 				NULL, 0);
    ConfigOutput(10, "ValidPRNs",			CHAR_t, 	VALUE,			CUDA_DEVICE, VECTORLENGTH_ANY, 	NULL, 0);
    ConfigOutput(11, "cpReference",			INT_t, 		VALUE,			CUDA_DEVICE, VECTORLENGTH_ANY, 	NULL, 0);
    ConfigOutput(12, "cpElapsedStart",		INT_t, 		VALUE,			CUDA_DEVICE, VECTORLENGTH_ANY, 	NULL, 0);
    ConfigOutput(13, "cpElapsedEnd",		INT_t, 		VALUE,			CUDA_DEVICE, VECTORLENGTH_ANY, 	NULL, 0);
    ConfigOutput(14, "ENU2ECEFMat",			DOUBLE_t,	VALUE,			CUDA_DEVICE, 9,  NULL, 0);
    ConfigOutput(15, "SatStatesOld", 			DOUBLE_t,	STATE,			CUDA_DEVICE, VECTORLENGTH_ANY,	NULL, 0);
    ConfigOutput(16, "cpRef",				INT_t, 		VALUE,			CUDA_DEVICE, VECTORLENGTH_ANY, 	NULL, 0);
    ConfigOutput(17, "cpRefTOW",		INT_t, 		VALUE,			CUDA_DEVICE, VECTORLENGTH_ANY, 	NULL, 0);
}


dsp::cuChanMgr::~cuChanMgr(){
	if (Started) Stop();
	delete [] ephSetArr;
    delete [] inputs;
    delete [] outputs;
    delete [] expectedInputs;
}


int
dsp::cuChanMgr::Start(void* cuFlowStream) {

	// Check module status and report accordingly
	if (Started) {
        std::clog << "[" << ModuleName << "] Start: Already Started." << std::endl;
        return 0;
    }
    std::clog << "[" << ModuleName << "] Starting ... " << std::flush;


	// Set the CUDA Stream for the GPU operations
    cuStream = (cudaStream_t*)cuFlowStream;



    // Get inputs
    // Determine how many channels are being tracked
    numChan = inputs[1]->VectorLength;

    // Get the pointer to the state (on device)
    xk1k1_d 	= (double*)inputs[10]->Data;
    xkk1_d		= (double*)inputs[12]->Data;

    // Get manifold info
    timeGrid_d 	= (double*)inputs[13]->Data;
    timeGridDim = inputs[13]->VectorLength;

    // Get the starting time
    rxTime	 	= *((double*)(inputs[9]->Data));

    // Length of one sample set (should be 1ms)
    T = *((double*)(inputs[11]->Data));
    // Round to the nearest 1us to ensure there is no trailing garbage precision
    T = round(T*1.0e6)/1.0e6;

    // Copy processing parameters to device constant memory
    cuCheckMSt(cudaMemcpyToSymbol(CHM_LCA_d, &numCA, sizeof(int), 0, cudaMemcpyHostToDevice));
    cuCheckMSt(cudaMemcpyToSymbol(CHM_T_d, &T, sizeof(double), 0, cudaMemcpyHostToDevice));



    // Get the input parameters copied into local device memory
    int size = sizeof(uint8_t)*CONST_PRN_MAX;
    cuCheckMSt(cudaMalloc((void**)&PRNs_d, size));
    cuCheckMSt(cudaMemcpyAsync(PRNs_d, (uint8_t*)inputs[1]->Data, size, cudaMemcpyHostToDevice, *cuStream));
    // Channel parameters and friends
    size = sizeof(double)*CONST_PRN_MAX;
    cuCheckMSt(cudaMalloc((void**)&codePhaseStart_d, size));
    cuCheckMSt(cudaMalloc((void**)&codePhaseEnd_d, size));
    cuCheckMSt(cudaMemcpyAsync(codePhaseEnd_d, (double*)inputs[2]->Data, size, cudaMemcpyHostToDevice, *cuStream));
    cuCheckMSt(cudaMalloc((void**)&carrierPhaseStart_d, size));
    cuCheckMSt(cudaMalloc((void**)&carrierPhaseEnd_d, size));
    cuCheckMSt(cudaMemcpyAsync(carrierPhaseEnd_d, (double*)inputs[3]->Data, size, cudaMemcpyHostToDevice, *cuStream));
    cuCheckMSt(cudaMalloc((void**)&codeFrequency_d, size));
    cuCheckMSt(cudaMemcpyAsync(codeFrequency_d, (double*)inputs[4]->Data, size, cudaMemcpyHostToDevice, *cuStream));
    cuCheckMSt(cudaMalloc((void**)&carrierFrequency_d, size));
    cuCheckMSt(cudaMemcpyAsync(carrierFrequency_d, (double*)inputs[5]->Data, size, cudaMemcpyHostToDevice, *cuStream));
    cuCheckMSt(cudaMalloc((void**)&txTime_d, size));
    size = sizeof(int)*CONST_PRN_MAX;
    cuCheckMSt(cudaMalloc((void**)&cpElapsedStart_d, size));
    cuCheckMSt(cudaMalloc((void**)&cpElapsedEnd_d, size));
    cuCheckMSt(cudaMemcpyAsync(cpElapsedEnd_d, (int*)inputs[6]->Data, size, cudaMemcpyHostToDevice, *cuStream));
    cuCheckMSt(cudaMalloc((void**)&cpReference_d, size));
    cuCheckMSt(cudaMemcpyAsync(cpReference_d, (int*)inputs[7]->Data, size, cudaMemcpyHostToDevice, *cuStream));
    cuCheckMSt(cudaMalloc((void**)&TOWcpReference_d, size));
    cuCheckMSt(cudaMemcpyAsync(TOWcpReference_d, (int*)inputs[8]->Data, size, cudaMemcpyHostToDevice, *cuStream));
    cuCheckMSt(cudaMalloc((void**)&dopplerSign_d, size));
    cuCheckMSt(cudaMemcpyAsync(dopplerSign_d, &dopplerSign, sizeof(int), cudaMemcpyHostToDevice, *cuStream));
    cuCheckMSt(cudaMalloc((void**)&ephToUse_d, size));
    // PVT state of each satellite as calculated by the channel params
    size = sizeof(dsp::utils::state_t<double>)*CONST_PRN_MAX;
    cuCheckMSt(cudaMalloc((void**)&satStates_d, size));
    // Satellite state things
    cuCheckMSt(cudaMalloc((void**)&batchSatStates_d, sizeof(dsp::utils::state_t<double>)*CONST_PRN_MAX*timeGridDim));
    cuCheckMSt(cudaMalloc((void**)&enu2ecefMat_d, sizeof(double)*9));





    // Load ephemeris
    ephSetVec 		= *((std::vector<dsp::utils::ephSet_t>*)(inputs[0]->Data));
    ephSetSize	 	= ephSetVec.size();

    // Convert eph's from vector to array so cudaMemcpy can easily copy to device
    ephSetArr = new dsp::utils::ephSet_t[ephSetSize];
    for (int i = 0; i < ephSetSize; i++) { ephSetArr[i].copyInto(ephSetVec[i]); }

    // Allocate space for every ephSet in navData
    size = sizeof(dsp::utils::ephSet_t)*ephSetSize;
    cuCheckMSt(cudaMalloc((void**)&ephSetPtr_d, size));
    cuCheckMSt(cudaMemcpyAsync(ephSetPtr_d, ephSetArr, size, cudaMemcpyHostToDevice));


    // Compute satellite states for the first computation
    // (launching kernels less than 32 threads can be more inefficient, even if only <10 channels are being tracked)
    CHM_ComputeSatStates<<<1, auxil::roundUpToNextPowerOfTwo(CONST_PRN_MAX), 0, *cuStream>>>
    		(ephSetPtr_d, ephSetSize, ephToUse_d,
    		 satStates_d,
    		 codePhaseEnd_d,
    		 cpElapsedEnd_d, cpReference_d, TOWcpReference_d,
    		 PRNs_d, numChan,
    		 txTime_d);

    // Propagate the channels for the next iteration
    // (launching kernels less than 32 threads can be more inefficient, even if only <10 channels are being tracked)
    CHM_TimeUpdateChannels<<<1, auxil::roundUpToNextPowerOfTwo(CONST_PRN_MAX), 0, *cuStream>>>
    		(ephSetPtr_d, ephSetSize, ephToUse_d,
    		 xk1k1_d, satStates_d,
    		 codePhaseStart_d, codePhaseEnd_d, codeFrequency_d,
    		 carrierPhaseStart_d, carrierPhaseEnd_d, carrierFrequency_d, dopplerSign_d,
    		 cpElapsedStart_d, cpElapsedEnd_d, cpReference_d, TOWcpReference_d,
    		 PRNs_d, numChan,
    		 rxTime, txTime_d, T);

    // Advance the rxTime to the end of the next sample set
    // (needs to happen after the propagate because MeasUpdate uses the rxTime of this iteration)
    rxTime += T;



    // Compute the rotated satellite positions
    prepBlockCount = floor((numChan*timeGridDim+1)/64)+1;
    // Shouldn't need to update sat states here since CHM_TimeUpdateChannels does it
    // Compute the satellite state from the ephemerides using the channel parameters
    // Rotate the satellites over the earth for every candidate grid point
    CHM_GridPrep<<<prepBlockCount, 64, 0, *cuStream>>> (rxTime, txTime_d, xkk1_d, satStates_d,
    		numChan, timeGrid_d, timeGridDim,
    		batchSatStates_d, enu2ecefMat_d);


    // Now that space is allocated, assign the position buffer and size
    outputs[0].Data = (void*)&rxTime;
    outputs[0].VectorLength = 1;
    outputs[1].Data = txTime_d;
    outputs[1].VectorLength = numChan;
    outputs[2].Data = codePhaseStart_d;
    outputs[2].VectorLength = numChan;
    outputs[3].Data = carrierPhaseStart_d;
    outputs[3].VectorLength = numChan;
    outputs[4].Data = codePhaseEnd_d;
	outputs[4].VectorLength = numChan;
	outputs[5].Data = carrierPhaseEnd_d;
	outputs[5].VectorLength = numChan;
    outputs[6].Data = codeFrequency_d;
	outputs[6].VectorLength = numChan;
	outputs[7].Data = carrierFrequency_d;
	outputs[7].VectorLength = numChan;
	outputs[8].Data = batchSatStates_d;
	outputs[8].VectorLength = numChan*timeGridDim;
	outputs[9].Data = dopplerSign_d;
	outputs[9].VectorLength = numChan;
	outputs[10].Data = PRNs_d;
	outputs[10].VectorLength = numChan;
	outputs[11].Data = cpReference_d;
	outputs[11].VectorLength = numChan;
	outputs[12].Data = cpElapsedStart_d;
	outputs[12].VectorLength = numChan;
	outputs[13].Data = cpElapsedEnd_d;
	outputs[13].VectorLength = numChan;
	outputs[14].Data = enu2ecefMat_d;
	outputs[14].VectorLength = 9;
	outputs[15].Data = satStates_d;
	outputs[15].VectorLength = numChan;
	outputs[16].Data = cpReference_d;
	outputs[16].VectorLength = numChan;
	outputs[17].Data = TOWcpReference_d;
	outputs[17].VectorLength = numChan;


    // Make sure all GPU tasks have completed before continuing
    cuCheckMSt(cudaStreamSynchronize(*cuStream));
    cuCheckMSt(cudaDeviceSynchronize());



    // Signifies that the next call to update() will be the first after start()
    Started = 1;

    std::clog << "Started." << std::endl;
    return 0;
}


int
dsp::cuChanMgr::Stop() {

    int ret = 0;
    if (Started == 0) {
        std::clog << "[" << ModuleName << "] Stop: Wasn't Started." << std::endl;
        return 0;
    }
    std::clog << "[" << ModuleName << "] Stopping ... " << std::flush;


    // Free device memory
    cuCheckMSp(cudaFree((void*)codePhaseStart_d));
    cuCheckMSp(cudaFree((void*)codePhaseEnd_d));
    cuCheckMSp(cudaFree((void*)carrierPhaseStart_d));
    cuCheckMSp(cudaFree((void*)carrierPhaseEnd_d));
    cuCheckMSp(cudaFree((void*)codeFrequency_d));
    cuCheckMSp(cudaFree((void*)carrierFrequency_d));
    cuCheckMSp(cudaFree((void*)cpElapsedStart_d));
    cuCheckMSp(cudaFree((void*)cpElapsedEnd_d));
    cuCheckMSp(cudaFree((void*)cpReference_d));
    cuCheckMSp(cudaFree((void*)TOWcpReference_d));
    cuCheckMSp(cudaFree((void*)dopplerSign_d));
    cuCheckMSp(cudaFree((void*)ephSetPtr_d));
    cuCheckMSp(cudaFree((void*)satStates_d));
    cuCheckMSp(cudaFree((void*)batchSatStates_d));
    cuCheckMSp(cudaFree((void*)enu2ecefMat_d));


    Started = 0;
    std::clog << "Stopped." << std::endl;

    return ret;
}


int
dsp::cuChanMgr::Update(void* cuFlowStream) {

    if (Started == 0){
        std::cerr << "[" << ModuleName
                  << "] Error: Update() Failed due to SatPos not initialized"
                  << std::endl;
        return -1;
    }


    // Propagate the channels for the next iteration
    // (launching kernels less than 32 threads can be more inefficient, even if only <10 channels are being tracked)
    CHM_PropagateChannels<<<1, auxil::roundUpToNextPowerOfTwo(CONST_PRN_MAX), 0, *cuStream>>>
    		(ephSetPtr_d, ephSetSize, ephToUse_d,
    		 xk1k1_d, satStates_d,
    		 codePhaseStart_d, codePhaseEnd_d, codeFrequency_d,
    		 carrierPhaseStart_d, carrierPhaseEnd_d, carrierFrequency_d, dopplerSign_d,
    		 cpElapsedStart_d, cpElapsedEnd_d, cpReference_d, TOWcpReference_d,
    		 PRNs_d, numChan,
    		 rxTime, txTime_d, T);


    // Advance the rxTime to the end of the next sample set
    // (needs to happen after the propagate because MeasUpdate uses the rxTime of this iteration)
    rxTime += T;



    // Compute the rotated satellite positions
    prepBlockCount = floor((numChan*timeGridDim+1)/64)+1;
    // Shouldn't need to run CHM_ComputeSatStates since CHM_PropagateChannels does this at the end
    // Determine the satellite positions from the ephemerides using the channel parameters
    // Rotate the satellite over the earth for each candidate grid point
    CHM_GridPrep<<<prepBlockCount, 64, 0, *cuStream>>> (rxTime, txTime_d, xkk1_d, satStates_d,
    		numChan, timeGrid_d, timeGridDim,
    		batchSatStates_d, enu2ecefMat_d);


    // Block on host until channel parameters have updated
    cuCheckMSt(cudaStreamSynchronize(*cuStream));


	return 0;
}


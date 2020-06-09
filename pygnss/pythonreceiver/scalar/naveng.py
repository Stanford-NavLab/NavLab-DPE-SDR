# -*- coding: utf-8 -*-



from ..libgnss import utils, satpos
from ..libgnss.constants import *

import numpy as np

def calculate_nav_soln(receiver, prn_list=None, mc=None, rxTime0=None, rxPos0=None, pOut=False):

    if prn_list is None:
        prn_list = sorted(receiver.channels.keys())

    if mc is None:
        mc = receiver._mcount

    numChannels = len(prn_list)

    codeIntDiff  = np.zeros(numChannels)
    codeFracDiff = np.zeros(numChannels)
    transmitTime = np.zeros(numChannels)

    sats_ECEF    = np.matrix(np.zeros((8,numChannels)))

    doppler      = np.zeros(numChannels)

    for prn_idx, prn in enumerate(prn_list):

        codeIntDiff[prn_idx]  = (receiver.channels[prn].cp[mc] \
                                 - receiver.channels[prn].ephemerides.timestamp['cp'])*T_CA
        codeFracDiff[prn_idx] = (receiver.channels[prn].rc[mc])/F_CA
        transmitTime[prn_idx] = receiver.channels[prn].ephemerides.timestamp['TOW'] \
                                 + codeIntDiff[prn_idx] + codeFracDiff[prn_idx]

        clkb, clkd = satpos.satellite_clock_correction(receiver.channels[prn].ephemerides,
                                                               transmitTime[prn_idx])
        sats_ECEF[:,prn_idx] = satpos.locate_satellite(receiver.channels[prn].ephemerides,
                                                               transmitTime[prn_idx]-clkb, clkb, clkd)

        doppler[prn_idx]      = (receiver.channels[prn].fi[mc])*receiver.rawfile.ds

    # Channels with a larger codeIntDiff + codeFracDiff are from satellites
    # which are closer in distance. Thus at a fixed received time instance,
    # the signals in those Channels correspond to a later transmitted time.

    # We need to remove satellite clock errors from the pseudorange
    # To reduce the number of for-loops, we determine that together with the
    # satellite positions in the Earth-Centered-Earth-Fixed frame (ECEF)

    if rxTime0 is None:
        # It takes about 68ms for the signal to travel from the closest
        # satellite to the receiver.
        rxTime = max(transmitTime)+0.068
        #print('Initializing rxTime as: %.3f'%(rxTime))
    else:
        rxTime = rxTime0

    # organize observations into pseudoranges and pseudorates
    pseudoranges = (C)*(rxTime - transmitTime) + (C)*(np.asarray(sats_ECEF)[3,:])
    pseudorates  = (-C/F_L1)*(doppler) + (C)*(np.asarray(sats_ECEF)[7,:])

    transmitTime = transmitTime - np.asarray(sats_ECEF)[3,:]

    # rotate satellite positions to Earth-Centered-Inertial (ECI)
    t_c = rxTime

    sats_ECI = np.matrix(np.zeros((8,numChannels)))

    for prn_idx, prn in enumerate(prn_list):
        sats_ECI[:,prn_idx] = utils.ECEF_to_ECI(sats_ECEF[:,prn_idx],t_gps=transmitTime[prn_idx], t_c=t_c)

    # perform navigation solution estimation through iterative least squares in ECI,
    # get actual GPS received time, rotate back to ECEF
    posvel_ECI  = perform_least_sqrs(sats_ECI, pseudoranges, pseudorates=pseudorates, rxPos0=rxPos0)
    rxTime_a    = rxTime - posvel_ECI[3,0]/C
    posvel_ECEF = utils.ECI_to_ECEF(posvel_ECI, t_gps=rxTime_a, t_c=t_c)

    # rotate all states to receiver's ECI
    t_c = rxTime_a
    posvel_ECI = utils.ECEF_to_ECI(posvel_ECEF, t_gps=rxTime_a, t_c=t_c)
    for prn_idx, prn in enumerate(prn_list):
        sats_ECI[:,prn_idx] = utils.ECEF_to_ECI(sats_ECEF[:,prn_idx], t_gps=transmitTime[prn_idx], t_c=t_c)

    if pOut == True:
        return rxTime_a, rxTime, posvel_ECEF, posvel_ECI, sats_ECI, pseudoranges, pseudorates

    return rxTime_a, rxTime, posvel_ECEF, posvel_ECI, sats_ECI

def get_satellite_positions(receiver, prn_list=None, mc=None, t_c=None):

    if prn_list is None:
        prn_list = sorted(receiver.channels.keys())

    if mc is None:
        mc = receiver._mcount

    numChannels = len(prn_list)

    codeIntDiff  = np.zeros(numChannels)
    codeFracDiff = np.zeros(numChannels)
    transmitTime = np.zeros(numChannels)

    sats_ECEF    = np.matrix(np.zeros((8,numChannels)))

    for prn_idx, prn in enumerate(prn_list):

        codeIntDiff[prn_idx]  = (receiver.channels[prn].cp[mc] \
                                 - receiver.channels[prn].ephemerides.timestamp['cp'])*T_CA
        codeFracDiff[prn_idx] = (receiver.channels[prn].rc[mc])/F_CA
        transmitTime[prn_idx] = receiver.channels[prn].ephemerides.timestamp['TOW'] \
                                 + codeIntDiff[prn_idx] + codeFracDiff[prn_idx]

        clkb, clkd = satpos.satellite_clock_correction(receiver.channels[prn].ephemerides,
                                                               transmitTime[prn_idx])
        sats_ECEF[:,prn_idx] = satpos.locate_satellite(receiver.channels[prn].ephemerides,
                                                               transmitTime[prn_idx]-clkb, clkb, clkd)
        transmitTime[prn_idx] = transmitTime[prn_idx] - clkb

    if t_c is not None:

        sats_ECI = np.matrix(np.zeros((8,numChannels)))

        for prn_idx, prn in enumerate(prn_list):
            sats_ECI[:,prn_idx] = utils.ECEF_to_ECI(sats_ECEF[:,prn_idx],
                                                    t_gps=transmitTime[prn_idx], t_c=t_c)

        return sats_ECI, transmitTime

    return sats_ECEF, transmitTime

def perform_least_sqrs(sats, pseudoranges, pseudorates=None, iterations=10, rxPos0=None):

    assert np.shape(sats)[0] == 8,\
         "Error: Expected satPos to be formatted as (8,N) np.matrix,\
         given satPos is formatted as:"+str(np.shape(satPos))
    satPos  = sats[0:3,:]
    satVel  = sats[4:7,:]
    numSats = np.shape(satPos)[1]

    assert np.shape(pseudoranges) == (numSats,), \
        "Error: Expected pseudoranges to be of shape ("+str(numSats)+",),\
         given pseudoranges has shape:"+str(np.shape(pseudoranges))

    if pseudorates is not None:
        assert np.shape(pseudorates) == (numSats,), \
            "Error: Expected pseudorates to be of shape ("+str(numSats)+",),\
             given pseudorates has shape:"+str(np.shape(pseudorates))

    # Initialize the helper variables.
    if rxPos0 is None:
        rxPos = np.matrix(np.zeros((4,1)))
    else:
        assert np.shape(rxPos0) == (4,1), \
             "Error: Expected rxPos0 to be of shape (4,1), \
             given rxPos0 has shape:"+str(np.shape(rxPos0))
        rxPos = rxPos0

    A = np.matrix(np.zeros((numSats,4)))
    A[:,3] = np.matrix(np.ones((numSats,1)))
    b = np.matrix(np.zeros((numSats,1)))

    # Iteratively get the receiver location.
    for i in range(iterations):

        # Gather information from each satellite to initialize the b vector and A matrix.
        for idx in range(numSats):

            # Fill in values to b vector.
            b[idx,0] = pseudoranges[idx] - (np.linalg.norm(satPos[:,idx]-rxPos[0:3,0]) + rxPos[3,0])

            # Fill in values to A matrix (linearized equation).
            A[idx,0:3] = (-(satPos[:,idx] - rxPos[0:3,0])/(np.linalg.norm(satPos[:,idx]-rxPos[0:3]))).T

        # If A is not full rank, we cannot proceed.
        if np.linalg.matrix_rank(A) != 4:
            print('Error: the linear least squares estimate of the receiver',
                  'position requires that the geometry A matrix be of rank 4. ',
                  'A is currently of rank: '+str(np.linalg.matrix_rank(A)))
            print('A matrix:'+str(A))

        # Run linear least squares to get the position update.
        x, residuals, rank, s = np.linalg.lstsq(A, b, rcond= None)

        # Now apply the position update.
        rxPos = rxPos + x
        #print('position_calculation',i,np.linalg.norm(x), residuals, rank)

        if np.linalg.norm(x) < 1.0e-7:
            break

    if np.linalg.norm(x) > 1.0e-5:
        print('Error: update is greater than 1.0e-5m after '+str(iterations)+' iterations.')

    # Initialize the helper variable.
    rxVel = np.matrix(np.zeros((4,1)))

    for idx in range(numSats):
        # Compute the updated los.
        los = ((satPos[:,idx] - rxPos[0:3,0])/(np.linalg.norm(satPos[:,idx]-rxPos[0:3]))).T

        # Fill in values to A matrix (unit Line-Of-Sight vectors)
        A[idx,0:3] = -los

        # Fill in values to the b vector
        b[idx] = pseudorates[idx] - los.dot(satVel[:,idx])

    # If A is not full rank, we cannot proceed and simply return zeros.
    if np.linalg.matrix_rank(A) != 4:
        print('Error: the linear least squares estimate of the receiver',
              'velocity requires that the geometry A matrix be of rank 4. ',
              'A is currently of rank: '+str(np.linalg.matrix_rank(A)))
        print('A matrix:'+str(A))

    # Run linear least squares to get the velocity.
    x, residuals, rank, s = np.linalg.lstsq(A, b, rcond = None)
    #print('velocity calculation',np.linalg.norm(x), residuals, rank)

    rxVel = x

    retMat = np.matrix(np.zeros((8,1)))
    retMat[0:4,0] = rxPos
    retMat[4:8,0] = rxVel
    return retMat

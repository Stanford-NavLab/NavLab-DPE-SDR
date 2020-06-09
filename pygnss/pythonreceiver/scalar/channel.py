# -*- coding: utf-8 -*-


from ..libgnss import lockdetector, snrmeter, utils, satpos
from ..libgnss.constants import *

import correlator, discriminator, loopfilter

import numpy as np
import math
import scipy.io as sio
import csv

# NOTE: cp_compl, cp_sign

# to log receiver synchronous measurements
_M_ARRAY_DATA_NAMES = \
    ['cp', 'rc', 'ri', 'fc', 'fi', 'iE', 'qE', 'iP', 'qP', 'iL', 'qL',
     'dc', 'di', 'efc', 'efi', 'dpc', 'dpi', 'dfc', 'dfi', 'fc_bias', 'fi_bias',
     'lock', 'lockval', 'snr',
     'cp_sign']

# bit_cp gives code period at the bit
# bit_sign gives the bit sign

_CHANNEL_ATTRIBUTE_NAMES    = ['prn','_cpcount']
_CORRELATOR_ATTRIBUTE_NAMES = ['p_a']

class Channel():
    """A receiver component for managing quantities related to a specific GNSS satellite in view.

    Channel interacts with channel-related objects through functions and attributes.
    Channel saves measurements as attribute arrays.
    """

    def __init__(self, prn, receiver, param_list=_M_ARRAY_DATA_NAMES):
        """
        Constructs a Channel object with specified prn and measurement length.
        """

        self.prn = prn
        self.receiver = receiver
        self.rawfile = self.receiver.rawfile

        self._measurement_logs_on = False
        self.init_measurement_logs(param_list=param_list)

        self._cpcount = 0 # internal codeperiod  counter
        self.cp[0] = 0

        # Initialize correlator.
        self.correlator = correlator.Correlator(self.prn, channel=self)

        # Initialize discriminators and loopfilters.
        self.cdiscriminator = discriminator.Discriminator(flavor='DLL', channel=self)
        self.idiscriminator = discriminator.Discriminator(flavor='PLL', channel=self)
        self.cloopfilter = loopfilter.LoopFilter(self.rawfile.T, Bnp=3.0, channel=self)
        self.iloopfilter = loopfilter.LoopFilter(self.rawfile.T, Bnp=40.0, channel=self)

        # Initialize lockdetector and snrmeter.
        self.lockdetector = lockdetector.LockDetector(N=20, k=1.5, lossthreshold=50, lockthreshold=240)
        self.snrmeter     = snrmeter.SignalNoiseMeter(N=20, T=self.rawfile.T)

    def set_params(self, rc=None, ri=None, fc=None, fi=None):
        """
        Code phase, carrier phase, code frequency, carrier doppler frequency.
        NOTE : loopfilter is reset, ok for scalar. incorrect for vector.
        """

        mc = self.receiver._mcount

        if ri is not None:
            self.ri[mc] = ri
        if fi is not None:
            self.fi[mc] = fi

        if rc is not None:
            self.rc[mc] = rc
        if fc is not None:
            self.fc[mc] = fc

    def set_scalar_params(self, rc=None, ri=None, fc=None, fi=None):
        """
        Code phase, carrier phase, code frequency, carrier doppler frequency.
        NOTE : loopfilter is reset, ok for scalar. incorrect for vector.
        """

        mc = self.receiver._mcount

        if ri is not None:
            self.ri[mc] = ri
        if fi is not None:
            self.fi[mc] = fi
            self.fi_bias[mc] = fi
            self.iloopfilter.reset()

        if rc is not None:
            self.rc[mc] = rc
        if fc is not None:
            self.fc[mc] = fc
            self.fc_bias[mc] = fc - F_CA - self.rawfile.fcaid*(self.fi_bias[mc])
            self.cloopfilter.reset()

    def scalar_correlation(self):

        mc = self.receiver._mcount

        # Correlation
        self.iE[mc], self.qE[mc], self.iP[mc], self.qP[mc], self.iL[mc], self.qL[mc], cp_compl, cp_sign \
        = self.correlator.scalar_correlate(self.rawfile, self.rc[mc], self.ri[mc], self.fc[mc], self.fi[mc])

        # Calculate lock and SNR
        self.lock[mc], self.lockval[mc] = self.lockdetector.update(self.iP[mc], self.qP[mc])
        self.snr[mc] = self.snrmeter.update(self.iP[mc], self.qP[mc])

        assert cp_compl in [0,1,2]

        for cp_compl_idx in range(cp_compl):
            self.cp_sign[self._cpcount] = cp_sign[cp_compl_idx]
            self._cpcount = self._cpcount + 1

        return

    def vector_correlation(self):

        mc = self.receiver._mcount

        # generate the in between shifts by shifting a little and performing correlation (FFT)

        self.code_corr, self.carr_fft, cp_compl = self.correlator.vector_correlate(self.rawfile,
                                                                                   self.rc[mc], self.ri[mc],
                                                                                   self.fc[mc], self.fi[mc],
                                                                                   self.cp[mc], self.ephemerides.timestamp['cp'])

        self._cpcount = self._cpcount + cp_compl

        return self.code_corr, self.carr_fft, cp_compl


    def vector_correlation_unfolded(self):

        mc = self.receiver._mcount

        # generate the in between shifts by shifting a little and performing correlation (FFT)

        self.code_corr, self.carr_fft, cp_compl = self.correlator.vector_correlate_unfolded(self.rawfile,
                                                                                   self.rc[mc], self.ri[mc],
                                                                                   self.fc[mc], self.fi[mc],
                                                                                   self.cp[mc], self.ephemerides.timestamp['cp'])

        self._cpcount = self._cpcount + cp_compl

        return self.code_corr, self.carr_fft, cp_compl




    def scalar_time_shift(self):
        """increment/step in time."""

        mc = self.receiver._mcount

        # Record number of completed code periods
        self._cpcount = self._cpcount \
            + int(np.floor((self.rawfile.S_skip*(self.fc[mc]/self.rawfile.fs)+self.rc[mc])/L_CA))

        # Increment phase
        self.rc[mc] = np.mod(self.rc[mc]+self.fc[mc]*self.rawfile.T_skip, L_CA)
        self.ri[mc] = np.mod(self.ri[mc]+self.fi[mc]*self.rawfile.T_skip, 1.0)

        return

    def scalar_time_update(self):
        """increment/step in time."""

        mc = self.receiver._mcount

        # Increment phase
        self.rc[mc+1] = np.mod(self.rc[mc]+self.fc[mc]*self.rawfile.T, L_CA)
        self.ri[mc+1] = np.mod(self.ri[mc]+self.fi[mc]*self.rawfile.T, 1.0)

        # Record number of completed code periods
        self.cp[mc+1] = self._cpcount

        # Frequency is the same
        self.fc[mc+1] = self.fc[mc]
        self.fi[mc+1] = self.fi[mc]
        self.fi_bias[mc+1] = self.fi_bias[mc]
        self.fc_bias[mc+1] = self.fc_bias[mc]

        return


    def scalar_time_update_adv(self, prn_idx):
        """increment/step in time."""

        mc = self.receiver._mcount

        # Frequency is the same
        self.fc[mc+1] = self.fc[mc]
        self.fi[mc+1] = self.fi[mc]
        self.fi_bias[mc+1] = self.fi_bias[mc]
        self.fc_bias[mc+1] = self.fc_bias[mc]

        # Increment carrier phase
        self.ri[mc+1] = np.mod(self.ri[mc]+self.fi[mc]*self.rawfile.T, 1.0)

        # Better (?) code phase update following Arthur's issue with results from code phase
        # Compute the expected code phase using the current channel params
        #cp_pred = self._cpcount = self._cpcount + int(
            #np.floor((self.rawfile.S_skip * (self.fc[mc] / self.rawfile.fs) + self.rc[mc]) / L_CA))

        cp_pred = self.cp[mc] + int(np.floor((self.rc[mc] + self.fc[mc] * self.rawfile.T)/ L_CA))
        rc_pred = np.mod(self.rc[mc] + self.fc[mc] * self.rawfile.T, L_CA)

        # Get the satellite position for the predicted txTime
        X_ECEF   = self.receiver.ekf.X_ECEF
        rxTime_a = self.receiver.rxTime_a
        rxTime = self.receiver.rxTime

        X_ECI  = utils.ECEF_to_ECI(X_ECEF,t_gps=rxTime_a,t_c=rxTime_a)

        transmitTime = self.ephemerides.timestamp['TOW'] \
                        + T_CA * (cp_pred - self.ephemerides.timestamp['cp']) \
                        + rc_pred / F_CA

        clkb, clkd = satpos.satellite_clock_correction(self.ephemerides, transmitTime)
        sats_ECEF = satpos.locate_satellite(self.ephemerides, transmitTime-clkb, clkb, clkd)
        transmitTime = transmitTime - clkb
        sats_ECI = utils.ECEF_to_ECI(sats_ECEF, t_gps=transmitTime, t_c=rxTime_a)


        bc_range = np.linalg.norm(sats_ECI[0:3] - X_ECI[0:3], axis=0)
        bc_pseudorange = bc_range + C * (X_ECI[3] / C - sats_ECI[3])
        bc_transmitTime = rxTime - bc_pseudorange / C
        bc_codeFracDiff = bc_transmitTime \
                          - self.ephemerides.timestamp['TOW'] \
                          - T_CA * (self.cp[mc] - self.ephemerides.timestamp['cp'])
        bc_rc = bc_codeFracDiff * F_CA

        # Compute the code phase for the current state
        self.cp[mc+1] = np.floor(bc_rc/L_CA) + self.cp[mc]
        self.rc[mc+1] = np.mod(bc_rc, L_CA)

        return

    def scalar_measurement_update(self):
        """Get measurements, correct channel parameters."""

        mc = self.receiver._mcount

        # Get carrier discriminations (errors)
        self.dpi[mc], self.dfi[mc] \
        = self.idiscriminator.update(self.iP[mc-1], self.qP[mc-1])

        # Get code discriminations (errors)
        self.dpc[mc], self.dfc[mc] \
        = self.cdiscriminator.update(self.iE[mc-1], self.qE[mc-1], self.iL[mc-1], self.qL[mc-1])

        # Filter, phase errors lumped in frequency
        self.di[mc] = self.iloopfilter.update(xp = self.dpi[mc], xf = self.dfi[mc])
        self.dc[mc] = self.cloopfilter.update(xp = self.dpc[mc], xf = self.dfc[mc])

        # Get code and carrier frequency errors, phase errors lumped in frequency
        self.efi[mc] =  (self.fi_bias[mc] + self.di[mc]) - self.fi[mc-1]
        self.efc[mc] = ((F_CA + self.fc_bias[mc] + self.dc[mc]) + self.rawfile.fcaid*(self.fi_bias[mc] + self.di[mc]))\
                       - self.fc[mc-1]

        # Correct carrier and code frequency
        self.fi[mc] = self.fi[mc-1]+ self.efi[mc]
        self.fc[mc] = self.fc[mc-1]+ self.efc[mc]

        return

    def init_measurement_logs(self, param_list):
        """Initialize temporary measurement logging arrays as attributes.

        Array values are initialized to np.nan.
        """

        assert not(self._measurement_logs_on), 'Please delete current measurement logs.'

        for name in param_list:
            setattr(self, name, np.ones(self.receiver.mcount_max)*np.nan)
        self._measurement_logs_on = True

    def delete_measurement_logs(self, param_list):
        """Delete current temporary measurement logging arrays."""

        assert self._measurement_logs_on, 'No measurement logs to be deleted.'

        for name in param_list:
            delattr(self, name)

        self._measurement_logs_on = False

    def save_measurement_logs(self, filename, csvfname, param_list=_M_ARRAY_DATA_NAMES):
        """
        Save measurement logs into a .mat file.
        """

        assert self._measurement_logs_on, 'No measurement logs to be saved.'

        save_dict = {}
        for name in param_list:
            save_dict['channel_array_'+name] = getattr(self,name)
        for name in _CHANNEL_ATTRIBUTE_NAMES:
            save_dict['channel_'+name] = getattr(self,name)
        for name in _CORRELATOR_ATTRIBUTE_NAMES:
            save_dict['correlator_'+name] = getattr(self.correlator,name)
        sio.savemat(filename,save_dict)

        with open(csvfname, 'wb') as csv_file:
            writer = csv.writer(csv_file)
            for key, value in save_dict.items():
                writer.writerow([key, value])

    def load_measurement_logs(self,filename, param_list=_M_ARRAY_DATA_NAMES):
        """
        Load measurement logs from a .mat file.

        note: scalar.LoopFilter history is not loaded, tracking results will deviate
              slightly at point of load. scalar.LoopFilter converges pretty quickly.
        """

        if self._measurement_logs_on:
            print('Warning: Overwriting existing measurement logs.')

        load_dict = sio.loadmat(filename)

        for name in param_list:
            if hasattr(self, name):
                tmpself = getattr(self,name)
                tmploadlen = len(load_dict['channel_array_'+name][0,:])
                tmpself[0:tmploadlen] = load_dict['channel_array_'+name][0,:]
                # loads from the front
            else:
                setattr(self,name,load_dict['channel_array_'+name][0,:])
        for name in _CHANNEL_ATTRIBUTE_NAMES:
            setattr(self,name,load_dict['channel_'+name][0,0])
        for name in _CORRELATOR_ATTRIBUTE_NAMES:
            setattr(self.correlator,name,load_dict['correlator_'+name][0,0])

    def get_Nms_correlation(self,ms,N):
        cp_idc  = np.mod(self.cp[(ms-N):(ms)]-self.ephemerides.timestamp['cp'],20)
        bd_idc, = np.where(np.diff(cp_idc)<0) #index of the last element of the previous navbit

        numBoundaries = np.shape(bd_idc)[0]
        assert(numBoundaries in [0,1,2]),"Error: Number of navigation bit boundaries found:"+str(numBoundaries)

        iE = self.iE[(ms-N):(ms)]
        iP = self.iP[(ms-N):(ms)]
        iL = self.iL[(ms-N):(ms)]
        qE = self.qE[(ms-N):(ms)]
        qP = self.qP[(ms-N):(ms)]
        qL = self.qL[(ms-N):(ms)]

        if numBoundaries == 0:
            return iE, iP, iL, qE, qP, qL

        combined = (iE+iP+iL)+1j*(qE+qP+qL)

        if numBoundaries == 0:
            return iE, iP, iL, qE, qP, qL

        combined = (iE+iP+iL)+1j*(qE+qP+qL)

        if numBoundaries == 1:

            # Segment the correlations into two parts
            idc_1 = range(0,(bd_idc[0]+1))
            idc_2 = range((bd_idc[0]+1),N)

            combined_1 = np.sum(combined[idc_1])
            combined_2 = np.sum(combined[idc_2])

            pos = np.abs(combined_1 + combined_2)
            neg = np.abs(combined_1 - combined_2)

            if neg > pos:
                iE[idc_2] = -iE[idc_2]
                iP[idc_2] = -iP[idc_2]
                iL[idc_2] = -iL[idc_2]
                qE[idc_2] = -qE[idc_2]
                qP[idc_2] = -qP[idc_2]
                qL[idc_2] = -qL[idc_2]


        else:
            # Segment the correlations into three parts
            idc_1 = range(0,(bd_idc[0]+1))
            idc_2 = range((bd_idc[0]+1),(bd_idc[1]+1))
            idc_3 = range((bd_idc[1]+1),N)

            combined_1 = np.sum(combined[idc_1])
            combined_2 = np.sum(combined[idc_2])
            combined_3 = np.sum(combined[idc_3])

            pos = np.abs(combined_1 + combined_2)
            neg = np.abs(combined_1 - combined_2)

            if neg > pos:
                iE[idc_2] = -iE[idc_2]
                iP[idc_2] = -iP[idc_2]
                iL[idc_2] = -iL[idc_2]
                qE[idc_2] = -qE[idc_2]
                qP[idc_2] = -qP[idc_2]
                qL[idc_2] = -qL[idc_2]
                combined_2 = -combined_2

            pos = np.abs(combined_1 + combined_2 + combined_3)
            neg = np.abs(combined_1 + combined_2 - combined_3)

            if neg > pos:
                iE[idc_3] = -iE[idc_3]
                iP[idc_3] = -iP[idc_3]
                iL[idc_3] = -iL[idc_3]
                qE[idc_3] = -qE[idc_3]
                qP[idc_3] = -qP[idc_3]
                qL[idc_3] = -qL[idc_3]

        return iE, iP, iL, qE, qP, qL


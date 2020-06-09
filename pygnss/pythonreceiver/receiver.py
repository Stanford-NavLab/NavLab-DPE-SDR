# -*- coding: utf-8 -*-

from vector import ekf
from scalar import channel, naveng, discriminator
from libgnss import dataparser, utils
from libgnss.constants import *

import os.path
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import csv

_RECEIVER_ATTRIBUTE_NAMES    = ['_mcount', 'mcount_max']
_RAWFILE_ATTRIBUTE_NAMES     = ['T', 'T_big']

_M_ARRAY_DATA_NAMES = \
    ['m_samp', 'm_time']
    # updated this to 'm_samp' and 'm_time', check

    #removed channel_list, change to channels.keys()



# Ephemerides attributes
# Note: this has to be in receiver so the ephemerides can be iterated through
# in order of prn list when logging to the hand-off file

_CLOCK_NAMES = ['weeknumber', 'accuracy', 'health', 'T_GD', 't_oc', 'a_f2',
                'a_f1', 'a_f0']

_EPHEMERIS_NAMES = ['C_rs', 'delta_n', 'M_0', 'C_uc', 'e', 'C_us',
                    'sqrt_A', 't_oe', 'C_ic', 'OMEGA_0', 'C_is', 'i_0',
                    'C_rc', 'omega', 'OMEGADOT', 'IDOT']

_NONSTANDARD_NAMES = ['total','complete','IODE','IODC']

_NAMES = _CLOCK_NAMES + _EPHEMERIS_NAMES + _NONSTANDARD_NAMES
_TOTAL = len(_NAMES)


class Receiver():
    """A PyGNSS receiver."""

    def __init__(self, rawfile, mcount_max=100000, channels = None):
        """Constructs a Receiver object that has a raw signal file.

        Args:
            rawfile (libgnss.rawfile.RawFile): A libgnss.rawfile.RawFile object.
            mcount_max (int): Maximum number of measurements to be performed.
                              Each measurement is associated with a certain time point/boundary.
            channels (dict): Dictionary of scalar.channel.Channel objects.
        """

        # Initialize link to rawfile
        self.rawfile = rawfile

        # Initialize measurement logs m_samp and m_time
        self.mcount_max = mcount_max
        self.init_measurement_logs(mcount_max)

        # Initialize dictionary of channels
        self.channels = channels if channels is not None else {}

        # Hold a reference mcount value
        self._ref_mcount = -1
        self._ref_bytes_read = -1

        return

    def init_measurement_logs(self, mcount_max):
        self._mcount = 0
        self.m_samp = np.ones(mcount_max)*np.nan
        self.m_time = np.ones(mcount_max)*np.nan
        return

    def add_channels(self, prn_list):

        for prn in prn_list:
            if prn in self.channels:
                print('Warning: Attempting to add PRN %d, overwriting.\n'%(prn))
            self.channels[prn] = channel.Channel(prn, self)

        return

    def del_channels(self, prn_list):

        for prn in prn_list:
            if prn in self.channels:
                del self.channels[prn]
            else:
                print('Warning: Attempting to remove nonexistent PRN %d.\n'%(prn))

        return

    def store_ref_mcount(self):
        self._ref_mcount = self._mcount
        self._ref_bytes_read = self.rawfile.rawfile.tell()
        return

    def init_dp(self, rxTime0=None, rxPos0=None, navguess = True, pOut=False):

        mc       = self._mcount
        prn_list = sorted(self.channels.keys())

        rxTime_a, rxTime, X_ECEF, X_ECI, sats_ECI = \
        naveng.calculate_nav_soln(self, mc=mc, prn_list=prn_list, rxTime0=rxTime0, rxPos0=rxPos0, pOut=pOut)

        self.ekf = ekf.ExtendedKalmanFilter(X_ECEF, Sigma=None, T=self.rawfile.T_big)
        self.rxTime = rxTime
        self.rxTime_a = rxTime_a

        # Debug -- remove the comment later
        #self.dp_measurement_update_channels()

        if navguess:
            self.navguess = NavigationGuesses()


    def init_handoff(self, rxTime0=None, rxPos0=None, navguess = True, pOut=False):
        mc = self._ref_mcount
        prn_list = sorted(self.channels.keys())

        rxTime_a, rxTime, X_ECEF, X_ECI, sats_ECI = \
        naveng.calculate_nav_soln(self, mc=mc, prn_list=prn_list, rxTime0=rxTime0, rxPos0=rxPos0, pOut=pOut)

        return rxTime, rxTime_a, X_ECEF

    def load_cudarecv_handoff(self, handoffFilename):
        with open (handoffFilename) as handoffFile:
            handoffReader = csv.reader(handoffFile, delimiter=',')

            for row in handoffReader:
                if   row[0] == 'rxTime':
                    self.rxTime = float(row[1])

                #elif row[0] == 'rxTime_a':
                #    self.rxTime_a = float(row[1])

                elif row[0] == 'X_ECEF':
                    X_ECEF = []
                    iterRow = iter(row)
                    next(iterRow)
                    for val in iterRow:
                        X_ECEF.append(float(val))
                    self.ekf.X_ECEF = np.matrix(X_ECEF).T
                    self.rxTime_a = self.rxTime - (self.ekf.X_ECEF[3]/C)

                elif row[0] == 'bytes_read':
                    self.rawfile.rawfile.seek(int(row[1]), 0)

                elif row[0] == 'rc':
                    for idx, prn in enumerate(self.channels):
                        self.channels[prn].rc[self._mcount] = float(row[idx+1])

                elif row[0] == 'ri':
                    for idx, prn in enumerate(self.channels):
                        self.channels[prn].ri[self._mcount] = float(row[idx+1])

                elif row[0] == 'fc':
                    for idx, prn in enumerate(self.channels):
                        self.channels[prn].fc[self._mcount] = float(row[idx+1])

                elif row[0] == 'fi':
                    for idx, prn in enumerate(self.channels):
                        self.channels[prn].fi[self._mcount] = float(row[idx+1])

                elif row[0] == 'cp':
                    for idx, prn in enumerate(self.channels):
                        self.channels[prn].cp[self._mcount] = float(row[idx+1])

                elif row[0] == 'cp_timestamp':
                    for idx, prn in enumerate(self.channels):
                        self.channels[prn].ephemerides.timestamp['cp'] = float(row[idx+1])

                elif row[0] == 'TOW':
                    for idx, prn in enumerate(self.channels):
                        self.channels[prn].ephemerides.timestamp['TOW'] = float(row[idx+1])


    def perturb_init_ENU(self, offset, refState):

        X_shift = utils.ENU_to_ECEF(refState, offset)
        self.ekf.X_ECEF[0:3] = X_shift

        return

    def perturb_init_ECEF(self, offset, refState):

        self.ekf.X_ECEF[0:3] = refState[0:3] + offset

        return



    def solve_scalar(self, rxTime0=None, rxPos0=None, navguess = True, pOut=False):
        mc = self._mcount
        prn_list = sorted(self.channels.keys())

        rxTime_a, rxTime, X_ECEF, X_ECI, sats_ECI = \
        naveng.calculate_nav_soln(self, mc=mc, prn_list=prn_list, rxTime0=rxTime0, rxPos0=rxPos0, pOut=pOut)

        return rxTime, rxTime_a, X_ECEF

    def dp_track(self, mtrack=3500):

        for i in range(mtrack):
            self.rawfile.seek_rawfile(self.rawfile.S_skip)
            self.rawfile.update_rawsnippet()
            self.m_samp[self._mcount] = self.rawfile.rawfile_samp
            self.m_time[self._mcount] = self.rawfile.rawfile_time

            self.dp_time_update_state()
            #self.dp_time_update_channels()
            self.dp_time_update_channels_unfolded()

            # 4. increment measurement count (change time point/boundary)
            self._mcount = self._mcount + 1

            #e = self.dp_measurement_estimation()
            e = self.dp_measurement_estimation_unfolded()
            self.dp_measurement_update_state(e)
            self.dp_measurement_update_channels()

        return

    def dp_time_update_state(self):

        self.ekf._time_update()
        self.rxTime   = self.rxTime+self.rawfile.T_big
        self.rxTime_a = self.rxTime-self.ekf.X_ECEF[3,0]/C

        return

    def dp_time_update_channels(self):

        for prn in self.channels:
            self.channels[prn].scalar_time_shift()
            self.channels[prn].vector_correlation()
            self.channels[prn].scalar_time_update()

        return

    def dp_time_update_channels_unfolded(self):

        for prn_idx, prn in enumerate(self.channels):
            self.channels[prn].scalar_time_shift()
            self.channels[prn].vector_correlation_unfolded()
            self.channels[prn].scalar_time_update_adv(prn_idx)

        return


    def dp_measurement_estimation(self,L_power=1,gXk_grid = None):

        mc       = self._mcount
        prn_list = sorted(self.channels.keys())

        X_ECEF   = self.ekf.X_ECEF
        rxTime_a = self.rxTime_a
        rxTime   = self.rxTime

        X_ECI  = utils.ECEF_to_ECI(X_ECEF,t_gps=rxTime_a,t_c=rxTime_a)
        sats_ECI, transmitTime = naveng.get_satellite_positions(self, t_c=rxTime_a)

        if gXk_grid is None: #single receiver mode
            gX_ECEF_fixed_vel, gX_ECI_fixed_vel, gX_ECEF_fixed_pos, gX_ECI_fixed_pos, bc_pos_corr, bc_vel_fft\
                = self.navguess.get_nav_guesses(X_ECEF, rxTime_a)
        else: # multi receiver mode
            gX_ECEF_fixed_vel, gX_ECEF_fixed_pos = gXk_grid # unpacking
            gX_ECI_fixed_vel = utils.ECEF_to_ECI(gX_ECEF_fixed_vel,t_gps=self.rxTime_a,t_c=self.rxTime_a)
            gX_ECI_fixed_pos = utils.ECEF_to_ECI(gX_ECEF_fixed_pos,t_gps=self.rxTime_a,t_c=self.rxTime_a)
            bc_pos_corr      = np.zeros(np.size(gX_ECEF_fixed_vel,1))
            bc_vel_fft       = np.zeros(np.size(gX_ECEF_fixed_pos,1))

        for prn_idx, prn in enumerate(prn_list):

            bc_los          = (sats_ECI[0:3,prn_idx]-X_ECI[0:3]).T/\
                              (np.linalg.norm((sats_ECI[0:3,prn_idx]-X_ECI[0:3]), axis=0)[0])
            bc_rangerate    = gX_ECI_fixed_pos[4:7,:]-sats_ECI[4:7,prn_idx]
            bc_losrangerate = np.array(bc_los*bc_rangerate)[0,:]
            bc_pseudorate   = -bc_losrangerate + C*(np.array(gX_ECI_fixed_pos)[7,:]/C-sats_ECI[7,prn_idx])
            bc_doppler      = -F_L1/C*bc_pseudorate
            bc_fi           = bc_doppler/self.rawfile.ds

            bc_fi0      = bc_fi - self.channels[prn].fi[mc]
            bc_fi0_idx  = (self.rawfile.carr_fftpts/self.rawfile.fs)*(bc_fi0)+self.rawfile.carr_fftpts/2.0
            bc_fi0_cidx = np.floor(bc_fi0_idx+1).astype(np.int32)
            bc_fi0_fidx = np.floor(bc_fi0_idx).astype(np.int32)
            bc_fi0_fft  = (self.channels[prn].carr_fft[bc_fi0_cidx]*(bc_fi0_idx-bc_fi0_fidx)\
                                +self.channels[prn].carr_fft[bc_fi0_fidx]*(bc_fi0_cidx-bc_fi0_idx))
            bc_vel_fft  = bc_vel_fft + np.abs(bc_fi0_fft)**L_power


            bc_range        = np.linalg.norm(sats_ECI[0:3,prn_idx]-gX_ECI_fixed_vel[0:3,:], axis=0)
            bc_pseudorange  = bc_range + C*(np.array(gX_ECI_fixed_vel)[3,:]/C-sats_ECI[3,prn_idx])
            bc_transmitTime = rxTime - bc_pseudorange/C
            bc_codeFracDiff = bc_transmitTime \
                              - self.channels[prn].ephemerides.timestamp['TOW']\
                              - T_CA*(self.channels[prn].cp[mc]-self.channels[prn].ephemerides.timestamp['cp'])
            bc_rc           = bc_codeFracDiff*F_CA

            bc_rc0      = bc_rc - self.channels[prn].rc[mc]

            bc_rc0_idx  = (self.rawfile.fs/self.channels[prn].fc[mc])*(-bc_rc0)+self.rawfile.S/(2.0*self.rawfile.N)
            bc_rc0_cidx = np.floor(bc_rc0_idx+1).astype(np.int32)
            bc_rc0_fidx = np.floor(bc_rc0_idx).astype(np.int32)
            bc_rc0_corr = (self.channels[prn].code_corr[bc_rc0_cidx]*(bc_rc0_idx-bc_rc0_fidx)\
                                +self.channels[prn].code_corr[bc_rc0_fidx]*(bc_rc0_cidx-bc_rc0_idx))
            bc_pos_corr = bc_pos_corr + np.abs(bc_rc0_corr)**L_power

        if gXk_grid is not None:
            self.pos_corr = bc_pos_corr
            self.vel_fft  = bc_vel_fft
            return

        mean_pos = np.sum(np.multiply(gX_ECEF_fixed_vel[0:4],np.tile(bc_pos_corr,(4,1))),axis=1)/np.sum(bc_pos_corr)
        mean_vel = np.sum(np.multiply(gX_ECEF_fixed_pos[4:8],np.tile(bc_vel_fft,(4,1))),axis=1)/np.sum(bc_vel_fft)

        return np.vstack((mean_pos,mean_vel))-X_ECEF




    def dp_measurement_estimation_unfolded(self,L_power=1,gXk_grid = None):

        mc       = self._mcount
        prn_list = sorted(self.channels.keys())

        X_ECEF   = self.ekf.X_ECEF
        rxTime_a = self.rxTime_a
        rxTime   = self.rxTime

        X_ECI  = utils.ECEF_to_ECI(X_ECEF,t_gps=rxTime_a,t_c=rxTime_a)
        sats_ECI, transmitTime = naveng.get_satellite_positions(self, t_c=rxTime_a)

        if gXk_grid is None: #single receiver mode
            gX_ECEF_fixed_vel, gX_ECI_fixed_vel, gX_ECEF_fixed_pos, gX_ECI_fixed_pos, bc_pos_corr, bc_vel_fft\
                = self.navguess.get_nav_guesses(X_ECEF, rxTime_a)
        else: # multi receiver mode
            gX_ECEF_fixed_vel, gX_ECEF_fixed_pos = gXk_grid # unpacking
            gX_ECI_fixed_vel = utils.ECEF_to_ECI(gX_ECEF_fixed_vel,t_gps=self.rxTime_a,t_c=self.rxTime_a)
            gX_ECI_fixed_pos = utils.ECEF_to_ECI(gX_ECEF_fixed_pos,t_gps=self.rxTime_a,t_c=self.rxTime_a)
            bc_pos_corr      = np.zeros(np.size(gX_ECEF_fixed_vel,1))
            bc_vel_fft       = np.zeros(np.size(gX_ECEF_fixed_pos,1))

        for prn_idx, prn in enumerate(prn_list):

            bc_los          = (sats_ECI[0:3,prn_idx]-X_ECI[0:3]).T/\
                              (np.linalg.norm((sats_ECI[0:3,prn_idx]-X_ECI[0:3]), axis=0)[0])
            bc_rangerate    = gX_ECI_fixed_pos[4:7,:]-sats_ECI[4:7,prn_idx]
            bc_losrangerate = np.array(bc_los*bc_rangerate)[0,:]
            bc_pseudorate   = -bc_losrangerate + C*(np.array(gX_ECI_fixed_pos)[7,:]/C-sats_ECI[7,prn_idx])
            bc_doppler      = -F_L1/C*bc_pseudorate
            bc_fi           = bc_doppler/self.rawfile.ds

            bc_fi0      = bc_fi - self.channels[prn].fi[mc]
            bc_fi0_idx  = (self.rawfile.carr_fftpts/self.rawfile.fs)*(bc_fi0)+self.rawfile.carr_fftpts/2.0
            #bc_fi0_idx = (self.rawfile.S / self.rawfile.fs) * (bc_fi0) + self.rawfile.S / 2.0
            #bc_fi0_idx = (131072 / self.rawfile.fs) * (bc_fi0) + 131072 / 2.0
            bc_fi0_cidx = np.floor(bc_fi0_idx+1).astype(np.int32)
            bc_fi0_fidx = np.floor(bc_fi0_idx).astype(np.int32)
            bc_fi0_fft  = (self.channels[prn].carr_fft[bc_fi0_cidx]*(bc_fi0_idx-bc_fi0_fidx)\
                                +self.channels[prn].carr_fft[bc_fi0_fidx]*(bc_fi0_cidx-bc_fi0_idx))
            bc_vel_fft  = bc_vel_fft + np.abs(bc_fi0_fft)**L_power


            bc_range        = np.linalg.norm(sats_ECI[0:3,prn_idx]-gX_ECI_fixed_vel[0:3,:], axis=0)
            bc_pseudorange  = bc_range + C*(np.array(gX_ECI_fixed_vel)[3,:]/C-sats_ECI[3,prn_idx])
            bc_transmitTime = rxTime - bc_pseudorange/C
            bc_codeFracDiff = bc_transmitTime \
                              - self.channels[prn].ephemerides.timestamp['TOW']\
                              - T_CA*(self.channels[prn].cp[mc]-self.channels[prn].ephemerides.timestamp['cp'])
            bc_rc           = bc_codeFracDiff*F_CA

            bc_rc0      = bc_rc - self.channels[prn].rc[mc]

            bc_rc0_idx  = (self.rawfile.fs/self.channels[prn].fc[mc])*(-bc_rc0)+self.rawfile.S/2.0
            bc_rc0_cidx = np.floor(bc_rc0_idx+1).astype(np.int32)
            bc_rc0_fidx = np.floor(bc_rc0_idx).astype(np.int32)
            bc_rc0_corr = (self.channels[prn].code_corr[bc_rc0_cidx]*(bc_rc0_idx-bc_rc0_fidx)\
                                +self.channels[prn].code_corr[bc_rc0_fidx]*(bc_rc0_cidx-bc_rc0_idx))
            bc_pos_corr = bc_pos_corr + np.abs(bc_rc0_corr)**L_power

        if gXk_grid is not None:
            self.pos_corr = bc_pos_corr
            self.vel_fft  = bc_vel_fft
            return

        #mean_pos = np.sum(np.multiply(gX_ECEF_fixed_vel[0:4],np.tile(bc_pos_corr,(4,1))),axis=1)/np.sum(bc_pos_corr)
        #mean_vel = np.sum(np.multiply(gX_ECEF_fixed_pos[4:8],np.tile(bc_vel_fft,(4,1))),axis=1)/np.sum(bc_vel_fft)

        mean_pos = gX_ECEF_fixed_vel[0:4, np.argmax(bc_pos_corr)]
        mean_vel = gX_ECEF_fixed_pos[4:8, np.argmax(bc_vel_fft)]

        temp = np.vstack((mean_pos,mean_vel))-X_ECEF
        return np.vstack((mean_pos,mean_vel))-X_ECEF




    def dp_measurement_update_state(self, e):

        self.ekf._measurement_update(e)
        self.rxTime_a = self.rxTime-self.ekf.X_ECEF[3,0]/C

        return


    # Ignore the first pass through here -- the result of dp_init() is overwritten
    def dp_measurement_update_channels(self):

        mc       = self._mcount
        prn_list = sorted(self.channels.keys())

        X_ECEF   = self.ekf.X_ECEF
        rxTime_a = self.rxTime_a
        rxTime   = self.rxTime

        X_ECI  = utils.ECEF_to_ECI(X_ECEF,t_gps=rxTime_a,t_c=rxTime_a)
        sats_ECI, transmitTime = naveng.get_satellite_positions(self, t_c=rxTime_a)

        for prn_idx, prn in enumerate(prn_list):

            bc_los          = (sats_ECI[0:3,prn_idx]-X_ECI[0:3]).T/\
                              (np.linalg.norm((sats_ECI[0:3,prn_idx]-X_ECI[0:3]), axis=0)[0])
            bc_rangerate    = X_ECI[4:7]-sats_ECI[4:7,prn_idx]
            bc_losrangerate = (bc_los*bc_rangerate)[0,0]
            bc_pseudorate   = -bc_losrangerate + C*(X_ECI[7,0]/C-sats_ECI[7,prn_idx])
            bc_doppler      = -F_L1/C*bc_pseudorate
            bc_fi           = bc_doppler/self.rawfile.ds
            #assert np.abs(bc_fi - self.channels[prn].fi[mc]) < 10, \
            #'prn: %d | bc_fi: %.3f, fi: %.3f'%(prn, bc_fi, self.channels[prn].fi[mc])
            self.channels[prn].fi[mc] = bc_fi

            bc_range         = np.linalg.norm(sats_ECI[0:3,prn_idx]-X_ECI[0:3], axis=0)[0]
            bc_pseudorange   = bc_range + C*(X_ECI[3,0]/C-sats_ECI[3,prn_idx])
            bc_transmitTime  = rxTime - bc_pseudorange/C
            bc_codeFracDiff  = bc_transmitTime \
                               - self.channels[prn].ephemerides.timestamp['TOW']\
                               - T_CA*(self.channels[prn].cp[mc]-self.channels[prn].ephemerides.timestamp['cp'])
            bc_rc            = bc_codeFracDiff*F_CA
            #assert np.abs(bc_rc - self.channels[prn].rc[mc]) < 1, \
            #'prn: %d | bc_rc: %.3f, rc: %.3f'%(prn, bc_rc, self.channels[prn].rc[mc])
            bc_fc            = F_CA + self.rawfile.fcaid*bc_fi + (bc_rc - self.channels[prn].rc[mc])/self.rawfile.T
            #assert np.abs(bc_fc - self.channels[prn].fc[mc]) < 20, \
            #'prn: %d | bc_fc: %.3f, fc: %.3f'%(prn, bc_fc, self.channels[prn].fc[mc])
            self.channels[prn].fc[mc] = bc_fc

        return

    def scalar_acquisition(self, prn_list, T=0.01, update=True):

        for prn in prn_list:
            if prn not in self.channels:
                print('Warning: Attempting to acquire nonexistent PRN %d, adding channel.\n'%(prn))
                self.add_channels([prn])

        original_T = self.rawfile.T
        original_T_big = self.rawfile.T_big

        acq_results = []

        self.rawfile.set_rawsnippet_settings(T=T, T_big=T)

        print('First %.3fs'%(T))
        self.rawfile.update_rawsnippet()
        self.m_samp[self._mcount] = self.rawfile.rawfile_samp
        self.m_time[self._mcount] = self.rawfile.rawfile_time

        for prn in prn_list:

            result_matrix, found, rc, ri, fc, fi, cppr, cppm = \
            self.channels[prn].correlator.search_signal(self.rawfile)

            acq_results.append([found, rc, ri, fc, fi, cppr, cppm])

            print('PRN: %d, Found: %s, Code: %.2f chips, Carrier: %.2f cycles, '
                  'Doppler: %.2f Hz, Cppr: %.2f, Cppm: %.2f'
                  %(prn, found, rc, ri, fi, cppr, cppm))

        print('')

        print('Second %.3fs'%(T))
        self.rawfile.update_rawsnippet()

        for prn_idx, prn in enumerate(prn_list):

            result_matrix, found, rc, ri, fc, fi, cppr, cppm = \
            self.channels[prn].correlator.search_signal(self.rawfile)

            if cppm > acq_results[prn_idx][-1]:

                rc = np.mod(rc-fc*T, L_CA)
                ri = np.mod(ri-fi*T, 1.0)

                if update:
                    self.channels[prn].set_scalar_params(rc=rc, ri=ri, fc=fc, fi=fi)

                print('PRN: %d, Found: %s, Code: %.2f chips, Carrier: %.2f cycles, '
                      'Doppler: %.2f Hz, Cppr: %.2f, Cppm: %.2f ***'
                      %(prn, found, rc, ri, fi, cppr, cppm))

            else:

                rc, ri, fc, fi = acq_results[prn_idx][1:5]

                if update:
                    self.channels[prn].set_scalar_params(rc=rc, ri=ri, fc=fc, fi=fi)

                print('PRN: %d, Found: %s, Code: %.2f chips, Carrier: %.2f cycles, '
                      'Doppler: %.2f Hz, Cppr: %.2f, Cppm: %.2f'
                      %(prn, found, rc, ri, fi, cppr, cppm))

        print('')

        self.rawfile.seek_rawfile(-2*self.rawfile.S)
        self.rawfile.set_rawsnippet_settings(T=original_T, T_big=original_T_big)

        return

    def scalar_track(self, mtrack = 36000):

        for i in range(mtrack):

            # 1. get rawsnippet
            self.rawfile.update_rawsnippet()
            self.m_samp[self._mcount] = self.rawfile.rawfile_samp
            self.m_time[self._mcount] = self.rawfile.rawfile_time

            # 2. perform correlation
            # 3. perform time update
            for prn in self.channels:
                self.channels[prn].scalar_correlation()
                self.channels[prn].scalar_time_update()

            # 4. increment measurement count (change time point/boundary)
            self._mcount = self._mcount + 1

            # 5. perform measurement update
            for prn in self.channels:
                self.channels[prn].scalar_measurement_update()


    def vt_init (self,numPrevSamples=20,N=20,rxPos0 = None, rxTime0 = None,pOut = False):#,N=20):
        mc       = self._mcount
        prn_list = sorted(self.channels.keys())

        rxTime_a, rxTime, self.X_ECEF, X_ECI, sats_ECI = \
        naveng.calculate_nav_soln(self, mc=mc, prn_list=prn_list, rxTime0=rxTime0, rxPos0=rxPos0, pOut=pOut)

        self.ekf = ekf.ExtendedKalmanFilter(self.X_ECEF, Sigma=None, T=self.rawfile.T_big,rxTime_a = rxTime_a,rxTime0 =rxTime)

        self.numPrevSamples = numPrevSamples
        self.N = N

        # generate container to store prev posvel_ECEF
        list_posvel_ECEF = []

        # initilize rxTime0 so that the clk bias has a standardized reference
        self.ekf.rxTime0  = np.round((self.ekf.rxTime_a - self.numPrevSamples * self.rawfile.T)*1000.0)/1000.0

        # perform position calculations for numPrevSamples x 20ms
        for ms in range(mc-20*20,mc,20):
            rxTime_a, rxTime, X_ECEF, X_ECI, sats_ECI = \
                    naveng.calculate_nav_soln(self, mc=ms, prn_list=prn_list, rxTime0=self.ekf.rxTime0, rxPos0 = None, pOut=False)

            list_posvel_ECEF.append(X_ECEF)
            self.ekf.rxTime0 += self.rawfile.T

        # calculate initSigma (ECEF)
        self.ekf.Sigma = np.matrix(np.cov(np.transpose(np.asarray(list_posvel_ECEF)[:,:,0])))
        print('initSigma', self.ekf.Sigma)

        for prn in self.channels.keys():
            self.channels[prn].idiscrimiantor = discriminator.Discriminator(flavor = 'FLL',channel = self.channels[prn])


    def vt_track (self,mtrack = 36000):
        for i in range(mtrack):
            # 1. get rawsnippet
            self.rawfile.update_rawsnippet()
            self.m_samp[self._mcount] = self.rawfile.rawfile_samp
            self.m_time[self._mcount] = self.rawfile.rawfile_time

            for prn in self.channels:
                self.channels[prn].scalar_correlation()

            self.vt_time_update()
            self._mcount += 1
            self.vt_measurement_update()


    def vt_measurement_update(self, channelList = None):
        ds = self.rawfile.ds
        N  = self.N

        if channelList is None:
            channelList = self.channels.keys()

        numChannels    = len(channelList)
        numPrevSamples = self.numPrevSamples

        mscount   = self._mcount
        ms_array  = np.arange(mscount-numPrevSamples,mscount)

        # CORRECTION: get error covariance W
        epc_var = np.matrix(np.zeros((numChannels,numPrevSamples))) # error in phase
        efi_var = np.matrix(np.zeros((numChannels,numPrevSamples))) # error in frequency

        # Looping through all the channels to get the error in phase and error in frequency (doppler shift)
        for idx, channelNr in enumerate(channelList):
            epc_var[idx, :] = channel[channelNr].dpc[ms_array-1]
            efi_var[idx, :] = channel[channelNr].dfi[ms_array-1]

        # Converting to SI units
        epc_var = epc_var*(-c/fc)    # adds to rc [chip * m/s * 1/(chips/s)]
        efi_var = efi_var*(-ds*c/L1) # adds to fi [cycles/s * m/s * 1/(cycles/s)]

        # Combining epc and efi into one vector and finding the variance (np.var)
        e_var = np.concatenate((epc_var,efi_var))
        e_var = np.var(e_var,axis=1)

        W = np.matrix(np.diag(np.asarray(e_var)[:,0]))

        # CORRECTION: get discriminations e

        epc = np.matrix(np.zeros((numChannels,1)))
        efi = np.matrix(np.zeros((numChannels,1)))

        # Actually finding the errors
        for channelNr in self.channel.keys():
            self.channel[channelNr].dpc[mscount-1], self.channel[channelNr].dfc[mscount-1] = \
                    dpc, dfc = self.channel[channelNr].cdiscriminator.update(self.iE[mc-1], self.qE[mc-1], self.iL[mc-1], self.qL[mc-1])

            self.channel[channelNr].dpi[mscount-1], self.channel[channelNr].dfi[mscount-1] = \
                    dpi, dfi = self.channel[channelNr].idiscriminator.update(self.iP[mc-2],self.qP[mc-2],self.iP[mc-1],self.qP[mc-1]) # carrier and

        # Input into error vector
        for idx, channelNr in enumerate(channelList):
            epc[idx] = self.channel[channelNr].dpc[mscount-1]
            efi[idx] = self.channel[channelNr].dfi[mscount-1]

        epc = epc*(-c/fc)    # adds to rc
        efi = efi*(-ds*c/L1) # adds to fi

        e  = np.concatenate((epc, efi))

        # UPDATE: get geometry matrix H (or Jacobian) (ECI) - basically one-step Newton Raphson

        sats_ECI, transmitTime = naveng.get_satellite_positions(self, prn_list = channelList, t_c= self.rxTime_a)
        X_ECI = libgnss.utils.ECEF_to_ECI(self.X_ECEF,t_gps=self.rxTime_a,t_c=self.rxTime_a)

        H = np.matrix(np.zeros((2*len(channelList),8)))

        # Populating the H matrix
        for idx, channelNr in enumerate(channelList):
            rxSat_range = np.linalg.norm(sats_ECI[0:3,idx]-X_ECI[0:3,0], axis=0)[0]
            rxSat_los   = ((sats_ECI[0:3,idx]-X_ECI[0:3,0])/rxSat_range).T
            H[idx,0:3]  = -rxSat_los; H[idx,3] = 1.0
            H[idx+len(channelList),4:7]  = -rxSat_los; H[idx+len(channelList),7] = 1.0

        K              = self.ekf._update_K(H,W)
        dX             = self.ekf._update_dX_ECI(e)
        X_ECEF, X_ECI  = self.ekf._correct_X_ECI() # Converts ECEF to ECI, updates X, then converts back to ECEF
        Sigma          = self.ekf._correct_Sigma(H)

        self.X_ECEF = X_ECEF

        return e, W, K, dX, X_ECEF, X_ECI, Sigma

    def vt_time_update(self):

        ds,fc,Tc,fcaid = self.rawfile.ds,F_CA,T_CA, ds * fc / F_L1
        N = self.N

        channel     = self.channel

        # UPDATE: fi, fc

        # Satellites in ECI again for calculations
        sats_ECI, transmitTime = naveng.get_satellite_positions(self, prn_list = channelList, t_c= self.rxTime_a)
        X_ECI = libgnss.utils.ECEF_to_ECI(self.X_ECEF,t_gps=self.rxTime_a,t_c=self.rxTime_a)
        X_ECI[3,0] = X_ECI[3,0]/c
        X_ECI[7,0] = X_ECI[7,0]/c

        # bc = "back computing"
        bc_doppler_all = np.zeros((numChannels))

        for inr, nr in enumerate(channelList):
            bc_los   = (sats_ECI[0:3,inr]-X_ECI[0:3,0]).T    # Line of sight
            bc_range = np.linalg.norm(bc_los, axis=1)[0];    # print(bc_range)  Range
            bc_los   = bc_los/bc_range                       # Line of sight unit vector
            bc_rangerate    = X_ECI[4:7,0]-sats_ECI[4:7,inr] # Range rate vector (time derivative of bc_los)
            bc_losrangerate = bc_los*bc_rangerate
            bc_pseudorate   = -bc_losrangerate[0,0] + c*(X_ECI[7,0]-sats_ECI[7,inr])
            bc_doppler      = -L1/c*bc_pseudorate
            bc_doppler_all[inr] = bc_doppler
            bc_fi           = bc_doppler/ds
            channel[nr]._fi = bc_fi

        for inr, nr in enumerate(channelList):
            bc_range = np.linalg.norm(sats_ECI[0:3,inr]-X_ECI[0:3,0], axis=0)[0]
            bc_pseudorange  = bc_range + c*(X_ECI[3,0]-sats_ECI[3,inr])
            bc_transmitTime = self.ekf.rxTime0 - bc_pseudorange/c
            channel[nr]._fc = fc + 1.0/(N*0.001)*((bc_transmitTime-(transmitTime[inr]+sats_ECI[3,inr]))*fc)\
                                 + fcaid*bc_doppler_all[inr]

        # PREDICT: new X, new Sigma
        self.X_ECEF   = self._predict_X()
        Q      = self.ekf._predict_Q()
        Sigma  = self.ekf._predict_Sigma()
        self.ekf.rxTime0  = np.round((self.ekf.rxTime0 + self.N*Tc)*1000.0)/1000.0
        self.ekf.rxTime_a = self.ekf.rxTime0 - self.X_ECEF[3,0]/c

        for i in range(N):
            for channelNr in self.channelList:
                channel[channelNr].update()

        return self.X_ECEF, Q, Sigma, self.ekf.rxTime0, self.rxTime_a


    def plot_qp_ip_correlations(self, prn_list = None, m_start = 0, m_end = None, ylim = None, xlim = None):

        if prn_list is None:
            prn_list = self.channels

        if m_end is None:
            m_end = self._mcount

        for prn in prn_list:
            plt.close('all')
            plt.plot(self.channels[prn].iP[m_start:m_end],self.channels[prn].qP[m_start:m_end],'b.')
            plt.title('Plot of quad-phase prompt correlation against in-phase prompt correlation'+
                      '\n PRN %d, from m_start = %d to m_end = %d'%(prn, m_start, m_end))
            plt.xlabel('in-phase prompt correlation')
            plt.ylabel('quad-phase prompt correlation')
            if (ylim is not None) and (xlim is not None):
                plt.ylim(ylim)
                plt.xlim(xlim)
            else:
                plt.axis('equal')
            plt.ticklabel_format(style='sci', axis='x', scilimits=(10,0))
            plt.ticklabel_format(style='sci', axis='y', scilimits=(10,0))
            plt.tight_layout()
            plt.show()

        return

    def plot_ip_correlations_time(self, prn_list = None, m_start = 0, m_end = None, ylim = None, xlim = None):

        if prn_list is None:
            prn_list = self.channels

        if m_end is None:
            m_end = self._mcount

        for prn in prn_list:
            plt.close('all')
            plt.plot(self.m_time[m_start:m_end],self.channels[prn].iP[m_start:m_end],'b.')
            plt.title('Plot of in-phase prompt correlation against time'+
                      '\n PRN %d'%(prn))
            plt.xlabel('time (s)')
            plt.ylabel('in-phase prompt correlation')
            if (ylim is not None) and (xlim is not None):
                plt.ylim(ylim)
                plt.xlim(xlim)
            #plt.ticklabel_format(style='sci', axis='x', scilimits=(10,0))
            plt.ticklabel_format(style='sci', axis='y', scilimits=(10,0))
            plt.tight_layout()
            plt.show()

        return

    def save_measurement_logs(self, dirname=None, subdir=''):
        """
        Save measurement logs into a .mat file.
        """

        if dirname is None:
            dirname = os.path.dirname(self.rawfile.abspath)

        save_dict = {}

        for name in _M_ARRAY_DATA_NAMES:
            save_dict['receiver_'+name] = getattr(self,name)
        for name in _RECEIVER_ATTRIBUTE_NAMES:
            save_dict['receiver_'+name] = getattr(self,name)
        for name in _RAWFILE_ATTRIBUTE_NAMES:
            save_dict['rawfile_'+name] = getattr(self.rawfile,name)

        save_dict['receiver_channels'] = sorted(self.channels.keys())
        if not os.path.exists(os.path.join(dirname,subdir)):
            os.makedirs(os.path.join(dirname,subdir))

        sio.savemat(os.path.join(dirname,subdir,'receiver.mat'),save_dict)


        for prn in self.channels:
            self.channels[prn].save_measurement_logs(os.path.join(dirname,subdir,'channel_%d.mat'%(prn)), os.path.join(dirname,subdir,'channel_%d.csv'%(prn)))

        return

    def save_scalar_handoff(self, prn_list, dirname=None, subdir=''):
        with open(os.path.join(dirname,subdir,'handoff_params.csv'), 'wb') as csv_file:
            writer = csv.writer(csv_file)

            # Record the receiver estimated position
            rxTime, rxTime_a, X_ECEF = self.init_handoff()

            writer.writerow(['rxTime', rxTime])
            writer.writerow(['rxTime_a', rxTime_a])
            x_guess = ['X_ECEF']

            for idx, val in enumerate(X_ECEF):
                x_guess.append(X_ECEF[idx, 0])

            writer.writerow(x_guess)

            # Record the file read position
            writer.writerow(['bytes_read', self._ref_bytes_read])

            prns = ['prn_list']
            rc = ['rc'] # code phase at mc
            ri = ['ri'] # carrier phase at mc
            fc = ['fc'] # code frequency at mc
            fi = ['fi'] # carrier frequency at mc
            cp = ['cp'] # code periods elpased at mc
            cp_timestamp = ['cp_timestamp'] # code periods elapsed at TOW time
            tow = ['TOW'] # reference TOW corresponding to cp_timestamp

            save_dict = {}


            for prn in self.channels:
                prns.append(prn)
                rc.append(self.channels[prn].rc[self._ref_mcount])
                ri.append(self.channels[prn].ri[self._ref_mcount])
                fc.append(self.channels[prn].fc[self._ref_mcount])
                fi.append(self.channels[prn].fi[self._ref_mcount])
                cp.append(self.channels[prn].cp[self._ref_mcount])

                if self.channels[prn].ephemerides is not None:
                    cp_timestamp.append(self.channels[prn].ephemerides.timestamp['cp'])
                    tow.append(self.channels[prn].ephemerides.timestamp['TOW'])
                else:
                    cp_timestamp.append(float('nan'))
                    tow.append(float('nan'))

                for name in _NAMES:
                    if self.channels[prn].ephemerides is not None:
                        if name in save_dict:
                            save_dict[name].append(getattr(self.channels[prn].ephemerides, name))
                        else:
                            save_dict[name] = [getattr(self.channels[prn].ephemerides, name)]
                    else:
                        save_dict[name].append(float('nan'))
                    # Iterate through the list to order params according to order of prn_list

            writer.writerow(prns)
            writer.writerow(rc)
            writer.writerow(ri)
            writer.writerow(fc)
            writer.writerow(fi)
            writer.writerow(cp)
            writer.writerow(cp_timestamp)
            writer.writerow(tow)

            for key, value in save_dict.items():
                temp = [key]
                temp.extend(value)
                writer.writerow(temp)

            print('Saved handoff params')
        return

    def load_measurement_logs(self, dirname=None, subdir='',prn_sel_list = None):
        """
        Load measurement logs from directory.

        note: scalar.LoopFilter history is not loaded, tracking results will deviate
              slightly at point of load. scalar.LoopFilter converges pretty quickly.
        """

        if dirname is None:
            dirname = os.path.dirname(self.rawfile.abspath)

        load_dict = sio.loadmat(os.path.join(dirname,subdir,'receiver.mat'))

        for name in _M_ARRAY_DATA_NAMES:
            setattr(self,name,load_dict['receiver_'+name][0,:])
        for name in _RECEIVER_ATTRIBUTE_NAMES:
            setattr(self,name,load_dict['receiver_'+name][0,0])
        for name in _RAWFILE_ATTRIBUTE_NAMES:
            setattr(self.rawfile,name,load_dict['rawfile_'+name][0,0])

        self.rawfile.set_rawsnippet_settings(T=self.rawfile.T, T_big=self.rawfile.T_big)

        if np.isnan(self.m_samp[self._mcount]):
            self.rawfile.seek_rawfile(self.m_samp[self._mcount-1]+self.rawfile.S,0)
        else:
            self.rawfile.seek_rawfile(self.m_samp[self._mcount],0)

        #if prn_list is None:
        prn_list = list(load_dict['receiver_channels'][0])
        # print 'prn_list',prn_list
        self.add_channels(prn_list)
        # print 'self.channels',self.channels
        del_list = []
        if prn_sel_list is not None:
            for prn in self.channels.keys():
                if prn not in prn_sel_list:
                    del_list +=[prn]
        print del_list
        self.del_channels(del_list)

        for prn in self.channels:
            self.channels[prn].load_measurement_logs(os.path.join(dirname,subdir,'channel_%d.mat'%(prn)))

        return

    def parse_ephemerides(self, prn_list=None, m_start=0, m_end=None):

        if prn_list is None:
            prn_list = self.channels

        if m_end is None:
            m_end = self._mcount-1

        for prn in prn_list:
            print('PRN %d'%(prn))
            dataparser.parse_ephemerides(self.channels[prn], m_start, m_end)

    def get_GDOP(self,X_ECEF = None,rxTime_a = None,prn_list = None,verbose = False):
        if prn_list is None:
            prn_list = self.channels

        if X_ECEF is None or rxTime_a is None:
            rxTime_a, rxTime, X_ECEF, X_ECI, sats_ECI = naveng.calculate_nav_soln(self,rxTime0 = rxTime_a,rxPos0 = X_ECEF)
        else:
            X_ECI  = utils.ECEF_to_ECI(X_ECEF,t_gps=rxTime_a,t_c=rxTime_a)
        sats_ECI, transmitTime = naveng.get_satellite_positions(self, prn_list = prn_list, t_c=rxTime_a)

        numChannels = len(prn_list)
        los = sats_ECI[0:3,:] - np.tile(X_ECI[0:3],(1,numChannels))
        los = (los/np.linalg.norm(los,axis=0)).T
        G    = np.hstack((-los,np.ones((numChannels,1))))
        G    = np.matrix(G)
        H    = np.linalg.inv(G.T*G)
        GDOP = np.sqrt(np.sum(np.diag(H)))
        if verbose:
            return GDOP,X_ECEF,rxTime_a,sats_ECI
        return GDOP

    def dp_get_GDOP(self,prn_list = None):
        return self.get_GDOP(X_ECEF = self.ekf.X_ECEF,rxTime_a = self.rxTime_a,prn_list = prn_list)

class NavigationGuesses():

    def __init__(self,gen = None):
        if gen:
            gen(self)
        else:
            #self.generate_evenly_spaced()
            self.generate_spread_grid()

    # predict the state X for the next timestamp
    def generate_evenly_spaced(self):

        dtmp  = np.linspace(-C/F_CA*2.0*0.6,C/F_CA*2.0*0.6,15)

        dZ = np.kron(dtmp,np.ones(15))
        dY = np.kron(dZ,np.ones(15))
        dX = np.kron(dY,np.ones(15))
        dY = np.tile(dY,15)
        dZ = np.tile(dZ,15*15)

        dXdot = dX/20.0
        dYdot = dY/20.0
        dZdot = dZ/20.0

        dT = np.linspace(-C/F_CA*2.0*0.6,C/F_CA*2.0*0.6,15)
        self.dT = np.tile(dT,15*15*15)

        dTdot = np.linspace(-C/F_L1*2.0*0.6,C/F_L1*2.0*0.6,15)
        self.dTdot = np.tile(dTdot,15*15*15)

        self.dX = np.mat(np.vstack((dX,dY,dZ)))
        self.dXdot = np.mat(np.vstack((dXdot,dYdot,dZdot)))

        self.N = len(self.dT)

        return
    
    def generate_spread_grid(self):

        dtmp = np.array([-22, -19, -16, -13, -10, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 10, 13, 16, 19, 22])

        dZ = np.kron(dtmp * 5,np.ones(25))
        dY = np.kron(dZ,np.ones(25))
        dX = np.kron(dY,np.ones(25))
        dY = np.tile(dY,25)
        dZ = np.tile(dZ,25*25)
        self.dX = np.mat(np.vstack((dX, dY, dZ)))

        dT = dtmp * 6
        self.dT = np.tile(dT, 25*25*25)


        dtmp = np.array([-12,-11,-10,-9,-8,-7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7,8,9,10,11,12])

        dZdot = np.kron(dtmp*5/10.0,np.ones(25))
        dYdot = np.kron(dZdot,np.ones(25))
        dXdot = np.kron(dYdot,np.ones(25))
        dYdot = np.tile(dYdot,25)
        dZdot = np.tile(dZdot,25*25)
        self.dXdot = np.mat(np.vstack((dXdot, dYdot, dZdot)))

        dTdot = dtmp * 0.25
        self.dTdot = np.tile(dTdot, 25*25*25)


        self.N = len(self.dT)
        self.Ndot = len(self.dTdot)

        return

    def get_nav_guesses(self, X_ECEF, rxTime_a, ECEF_only=False, scale = 1.0):

        gX_ECEF  = np.matrix(np.zeros((8,self.N)))

        X_ENU_trash, R = utils.ECEF_to_ENU(refState=X_ECEF[0:3], curState=X_ECEF[0:3])
        gX_ECEF[0:3,:] = utils.ENU_to_ECEF(refState=X_ECEF[0:3], diffState=self.dX    * scale, R_ECEF2ENU=R)
        gX_ECEF[3,:]   = X_ECEF[3,0]+self.dT
        gX_ECEF[4:7,:] = utils.ENU_to_ECEF(refState=X_ECEF[4:7], diffState=self.dXdot, R_ECEF2ENU=R)
        gX_ECEF[7,:]   = X_ECEF[7,0]+self.dTdot

        gX_ECEF_fixed_vel = np.matrix(gX_ECEF)
        gX_ECEF_fixed_vel[4:,:] = X_ECEF[4:,:]

        gX_ECEF_fixed_pos = np.matrix(gX_ECEF)
        gX_ECEF_fixed_pos[:4,:] = X_ECEF[:4,:]

        if ECEF_only:
            return gX_ECEF_fixed_vel, gX_ECEF_fixed_pos

        # State in ECI, fixed vel
        gX_ECI_fixed_vel = utils.ECEF_to_ECI(gX_ECEF_fixed_vel,t_gps=rxTime_a,t_c=rxTime_a)
        # State in ECI, fixed pos
        gX_ECI_fixed_pos = utils.ECEF_to_ECI(gX_ECEF_fixed_pos,t_gps=rxTime_a,t_c=rxTime_a)

        bc_pos_corr = np.zeros(self.N)
        bc_vel_fft  = np.zeros(self.N)

        return gX_ECEF_fixed_vel, gX_ECI_fixed_vel, gX_ECEF_fixed_pos, gX_ECI_fixed_pos, bc_pos_corr, bc_vel_fft

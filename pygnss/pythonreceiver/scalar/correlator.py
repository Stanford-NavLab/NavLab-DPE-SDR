# -*- coding: utf-8 -*-


from ..libgnss.constants import *

import numpy as np
import scipy.stats as sps
import math

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

B = 25;  dHz = 500; DOPPLER_SEARCH_MATRIX_NONCOHERENT = np.mat(np.arange((1-B)/2, (B-1)/2 + 1)).T*dHz
B = 125; dHz = 100; DOPPLER_SEARCH_MATRIX_COHERENT    = np.mat(np.arange((1-B)/2, (B-1)/2 + 1)).T*dHz

class Correlator():
    """A channel component for performing signal correlations with given raw signal snippet.

    Correlator searches for signal via brute force correlations.
    Correlator performs correlations against given raw signal snippet.
    """

    def __init__(self, prn, channel=None):
        """Constructs a Correlator object with specified prn.
        """

        self.prn = prn

        if channel is not None:
            self.channel = channel
            assert self.prn == self.channel.prn

        self.chips = self._make_L1_CAcode_chips()
        self.offset = 0.5

        self.p_a = 0

    def search_signal(self,
                      rawfile,
                      doppler_search_matrix = DOPPLER_SEARCH_MATRIX_COHERENT,
                      coherent = True,
                      plot_signal = False):
        """Calls both coarse_acquisition and fine_frequency_acquisition.
        """

        coarse_result_matrix, found, rc, fc, fi, cppr, cppm = self.coarse_acquisition(
                                                              rawfile, doppler_search_matrix, coherent, plot_signal)

        rc, ri, fc, fi = self.fine_frequency_acquisition(rawfile, rc, fc, fi, doppler_search_matrix)

        return coarse_result_matrix, found, rc, ri, fc, fi, cppr, cppm

    def coarse_acquisition(self, rawfile, doppler_search_matrix, coherent, plot_signal):
        """
        Creates 2D coarse_result_matrix,
           dimension 0 = doppler freq
           dimension 1 = code chips

        Returns coarse_result_matrix, found, rc, fc, fi, cppr, cppm
        """

        # Create 2D doppler wipeoff signal matrix
        doppler_wipeoff_signal_matrix = np.array(np.exp(-1j*(2*PI*doppler_search_matrix*rawfile.time_idc)))

        # Create code replica signal
        code_replica_signal = self.chips.take(np.mod(np.floor(rawfile.code_idc),L_CA).astype(np.int32))
        code_replica_signal_cfft = np.conjugate(np.fft.fft(code_replica_signal))

        # Initialize 2D coarse result matrix
        coarse_result_matrix = np.zeros(np.shape(doppler_wipeoff_signal_matrix),dtype=complex)

        # Perform brute-force correlations across dopplers
        for doppler_idx, doppler_signal in enumerate(doppler_wipeoff_signal_matrix):
            code_signal_fft = np.fft.fft(rawfile.rawsnippet*doppler_signal)
            coarse_result_matrix[doppler_idx,:] = np.fft.ifft(code_signal_fft*code_replica_signal_cfft)

        if rawfile.N != 1:
            tmp = np.reshape(coarse_result_matrix,(np.shape(coarse_result_matrix)[0],rawfile.N,rawfile.S/rawfile.N))
            if coherent :
                coarse_result_matrix = np.sum(tmp,axis=1)
            else:
                coarse_result_matrix = np.sum(np.abs(tmp),axis=1)

        coarse_result_matrix_abs = np.abs(coarse_result_matrix)

        # Search for maximum correlation, estimate code phase, carrier doppler frequency
        max_percode  = np.max(coarse_result_matrix_abs,axis=0)
        max_code_idx = max_percode.argmax()
        max_dopp_idx = coarse_result_matrix_abs[:,max_code_idx].argmax()
        rc = L_CA-rawfile.code_idc[max_code_idx]
        fi = doppler_search_matrix[max_dopp_idx,0]
        fc = F_CA + rawfile.fcaid*fi

        # Mask maximum correlation, estimate peak statistics
        peak = max_percode[max_code_idx]
        mask_S = int(np.ceil(rawfile.fs/F_CA))
        mask_idx = np.arange(-mask_S, mask_S+1) + max_code_idx

        max_percode[mask_idx] = 0
        cppr = peak / np.max(max_percode)
        cppm = peak / self._trim_mean(max_percode, 10)

        return coarse_result_matrix, cppm>2.0, rc, fc, fi, cppr, cppm

    def fine_frequency_acquisition(self, rawfile, rc, fc, fi, doppler_search_matrix):
        """
        Returns ri, fi
        """

        # Generate code replica signal
        # Update code_idc based on new fc, generate code replica signal
        # Doesn't seem to make a difference here

        code_idc = rawfile.time_idc*fc
        code_replica_signal = self.chips.take(np.mod(np.floor(code_idc+rc),L_CA).astype(np.int32))

        # Code wipeoff to get carrier only signal
        carr_signal = (rawfile.rawsnippet-np.mean(rawfile.rawsnippet)) * code_replica_signal

        # FFT to get carrier in frequency domain
        carr_signal_fft = np.fft.fftshift(np.fft.fft(carr_signal,rawfile.carr_fftpts))

        # Mask frequencies outside search range
        carr_signal_fft[rawfile.carr_fftidc<np.min(doppler_search_matrix)] = 0.0
        carr_signal_fft[rawfile.carr_fftidc>np.max(doppler_search_matrix)] = 0.0

        max_carr_idx = np.abs(carr_signal_fft).argmax()
        fine_frequency_result = carr_signal_fft[max_carr_idx]
        ri = np.angle(fine_frequency_result)/(2.0*PI)
        fi = rawfile.carr_fftidc[max_carr_idx]
        fc = F_CA + rawfile.fcaid*fi

        return rc, ri, fc, fi

    def scalar_correlate(self, rawfile, rc, ri, fc, fi):
        """Correlate received signal against replica signal.
        """

        # Carrier wipeoff
        baseband = rawfile.rawsnippet*np.exp(-1j*((2.0*PI*fi*rawfile.time_idc) + (2.0*PI*ri)))

        # Produce indices
        fidc = rawfile.time_idc*fc + rc

        eidc = np.mod(np.floor(fidc + self.offset),L_CA).astype(np.int32) # early indices
        pidc = np.mod(np.floor(fidc              ),L_CA).astype(np.int32) # prompt indices
        lidc = np.mod(np.floor(fidc - self.offset),L_CA).astype(np.int32) # late indices

        # Produce shifted replica signals.
        # numpy.take returns array values from given indices
        early  = self.chips.take(eidc)
        prompt = self.chips.take(pidc)
        late   = self.chips.take(lidc)

        # Determine boundaries in samples
        # (+1 is added to include the sample during array indexing)
        idxs1 = np.floor((L_CA-rc) * (rawfile.fs/fc)).astype(np.int32) + 1
        idxs2 = np.floor((2.0*L_CA-rc) * (rawfile.fs/fc)).astype(np.int32) + 1

        # The normal case: the first boundary is within this window, and the
        # second boundary is outside this window.
        if  idxs1 <= rawfile.S < idxs2:
            # The flow of this case is as follows:
            # First correlate part B and form signal   synchronous outputs.
            # Next  correlate part A and form receiver synchronous outputs.
            # |bbb|aaaaaa| |bbb|aaaaaa|

            # Segment the baseband vector into B and A parts.
            baseband_b = baseband[:idxs1]
            baseband_a = baseband[idxs1:]

            # Get the correlations for segments B and A.
            e_b = np.inner(baseband_b, early[:idxs1])  # early  part B
            p_b = np.inner(baseband_b, prompt[:idxs1]) # prompt part B
            l_b = np.inner(baseband_b, late[:idxs1])   # late   part B
            e_a = np.inner(baseband_a, early[idxs1:])  # early  part A
            p_a = np.inner(baseband_a, prompt[idxs1:]) # prompt part A
            l_a = np.inner(baseband_a, late[idxs1:])   # late   part A

            # Prepare the signal synchronous correlation outputs.
            p_s1 = self.p_a + p_b # prompt signal synchronous

            # Record the part A correlations for use in the next update.
            self.p_a = p_a

            # Prepare the receiver synchronous correlation outputs.
            pos = np.abs(e_b + p_b + l_b + e_a + p_a + l_a)
            neg = np.abs(e_b + p_b + l_b - e_a - p_a - l_a)

            e_r, p_r, l_r = 0, 0, 0

            if pos > neg: # there was no polarity change from B to A
                e_r = e_b + e_a # early  receiver synchronous
                p_r = p_b + p_a # prompt receiver synchronous
                l_r = l_b + l_a # late   receiver synchronous
            else: # there was a polarity change from B to A
                e_r = e_b - e_a # early  receiver synchronous
                p_r = p_b - p_a # prompt receiver synchronous
                l_r = l_b - l_a # late   receiver synchronous

            return e_r.real, e_r.imag, p_r.real, p_r.imag, l_r.real, l_r.imag, 1, -np.sign([p_s1.real])

        if  idxs1 < idxs2 <= rawfile.S:

            # The flow of this case is as follows:
            # First correlate part B and form signal synchronous outputs.
            # Next correlate an A & B pair and form signal synchronous outputs.
            # Finally correlate part A and form receiver synchronous outputs.

            # Segment the baseband vector into B, AB, and A parts.
            baseband_b = baseband[:idxs1]
            baseband_s = baseband[idxs1:idxs2]
            baseband_a = baseband[idxs2:]

            # Get the correlations for the first B segment.
            e_b = np.inner(baseband_b, early[:idxs1])  # early  part B
            p_b = np.inner(baseband_b, prompt[:idxs1]) # prompt part B
            l_b = np.inner(baseband_b, late[:idxs1])   # late   part B

            # Prepare the signal synchronous correlation outputs.
            p_s1 = self.p_a + p_b

            # Get the correlations for the complete A and B segment pair.
            e_s = np.inner(baseband_s, early[idxs1:idxs2])  # early  part S
            p_s = np.inner(baseband_s, prompt[idxs1:idxs2]) # prompt part S
            l_s = np.inner(baseband_s, late[idxs1:idxs2])   # late   part S

            p_s2 = p_s

            # Get the correlations for the final A segment.
            e_a = np.inner(baseband_a, early[idxs2:])  # early  part A
            p_a = np.inner(baseband_a, prompt[idxs2:]) # prompt part A
            l_a = np.inner(baseband_a, late[idxs2:])   # late   part A

            # Record the part A correlations for use in the next update.
            self.p_a = p_a

            # Prepare the receiver synchronous correlation outputs.
            pos = np.abs(e_b + p_b + l_b + e_s + p_s + l_s)
            neg = np.abs(e_b + p_b + l_b - e_s - p_s - l_s)

            e_r, p_r, l_r = 0, 0, 0

            if pos > neg: # there was no polarity change from B to S

                pos = np.abs(e_s + p_s + l_s + e_a + p_a + l_a)
                neg = np.abs(e_s + p_s + l_s - e_a - p_a - l_a)

                if pos > neg: # there was no polarity change from S to A
                    e_r = e_b + e_s + e_a
                    p_r = p_b + p_s + p_a
                    l_r = l_b + l_s + l_a
                else: # there was a polarity change from S to A
                    e_r = e_b + e_s - e_a
                    p_r = p_b + p_s - p_a
                    l_r = l_b + l_s - l_a

            else: # there was a polarity change from B to S (thus not S to A)

                e_r = e_b - e_s - e_a
                p_r = p_b - p_s - p_a
                l_r = l_b - l_s - l_a

            return e_r.real, e_r.imag, p_r.real, p_r.imag, l_r.real, l_r.imag, 2, -np.sign([p_s1.real,p_s2.real])

        if  rawfile.S < idxs1:
            # The flow of this case is as follows:
            # First correlate part B, actually a continuation of A from the
            # last correlation window.  Form the receiver synchronous outputs.

            # Get the correlations for the continuous B segment.
            e_b = np.inner(baseband, early)
            p_b = np.inner(baseband, prompt)
            l_b = np.inner(baseband, late)

            # Update the part A correlations for use in the next update.
            self.p_a = self.p_a + p_b

            return e_b.real, e_b.imag, p_b.real, p_b.imag, l_b.real, l_b.imag, 0, []

        else:
            print('EXTREME ERROR in scalar correlator')
            return

    def vector_correlate(self, rawfile, rc, ri, fc, fi, cp, cp_timestamp):

        cp_since_prev_bit = (cp-cp_timestamp)%20
        cp_to_next_bit    = 20-cp_since_prev_bit
        idx_next_bit      = np.floor((L_CA*cp_to_next_bit-rc)*(rawfile.fs/fc)).astype(np.int32) + 1

        cp_compl = np.floor((rawfile.S*(fc/rawfile.fs)+rc)/L_CA)

        if (idx_next_bit > 0) and (idx_next_bit < rawfile.S):

            raw_no_flip = np.array(rawfile.rawsnippet)
            raw_flip    = np.array(raw_no_flip)
            raw_flip[idx_next_bit:] = -raw_flip[idx_next_bit:]

            # #### convert to baseband with current fi and ri (velocities)

            doppler_wipeoff_signal = np.exp(-1j*((2.0*PI*fi*rawfile.time_idc) + (2.0*PI*ri)))

            baseband_no_flip = raw_no_flip*doppler_wipeoff_signal
            baseband_flip    = raw_flip*doppler_wipeoff_signal

            # #### already in baseband with current fi and ri (velocities)
            # - can do correlation (use shifted replica, will produce code phase error)
            # - can do code wipeoff and fft (will produce carrier freq error)

            # ##### correlation, save code_corr

            r = self.chips.take(np.mod(np.floor(rawfile.time_idc*fc + rc),L_CA).astype(np.int32))

            rcfft = np.conjugate(np.fft.fft(r))

            rfft_no_flip = np.fft.fft(baseband_no_flip)
            corr_no_flip = np.fft.ifft(rcfft*rfft_no_flip)
            corr_no_flip = np.sum(np.reshape(corr_no_flip,(rawfile.N,rawfile.S/rawfile.N)),axis=0)/rawfile.N

            rfft_flip = np.fft.fft(baseband_flip)
            corr_flip = np.fft.ifft(rcfft*rfft_flip)
            corr_flip = np.sum(np.reshape(corr_flip,(rawfile.N,rawfile.S/rawfile.N)),axis=0)/rawfile.N

            if np.abs(corr_flip[0])>np.abs(corr_no_flip[0]):
                raw = raw_flip
                code_corr = np.fft.fftshift(corr_flip)
            else:
                raw = raw_no_flip
                code_corr = np.fft.fftshift(corr_no_flip)

            # ##### FFT, save carr_fft, carr_fftbins

            baseband = (raw-np.mean(raw)) * r * doppler_wipeoff_signal
            carr_fft = np.fft.fftshift(np.fft.fft(baseband, rawfile.carr_fftpts))

        else:

            raw = np.array(rawfile.rawsnippet)
            doppler_wipeoff_signal = np.exp(-1j*((2.0*PI*fi*rawfile.time_idc) + (2.0*PI*ri)))

            baseband = raw*doppler_wipeoff_signal

            r = self.chips.take(np.mod(np.floor(rawfile.time_idc*fc + rc),L_CA).astype(np.int32))
            rcfft = np.conjugate(np.fft.fft(r))

            rfft = np.fft.fft(baseband)
            corr = np.fft.ifft(rcfft*rfft)
            corr = np.sum(np.reshape(corr,(rawfile.N,rawfile.S/rawfile.N)),axis=0)/rawfile.N

            code_corr    = np.fft.fftshift(corr)

            baseband = (raw-np.mean(raw)) * r * doppler_wipeoff_signal
            carr_fft = np.fft.fftshift(np.fft.fft(baseband, rawfile.carr_fftpts))

        while np.max(np.abs(carr_fft)) < np.max(np.abs(code_corr))/2:
            print('WARNING')

            max_code_idx = (np.abs(code_corr).argmax()-rawfile.S/(2.0*rawfile.N)).astype(np.int32)
            baseband = (raw-np.mean(raw)) * np.roll(r,max_code_idx) * doppler_wipeoff_signal
            carr_fft = np.fft.fftshift(np.fft.fft(baseband, rawfile.carr_fftpts))

        return code_corr, carr_fft, cp_compl




    def vector_correlate_unfolded(self, rawfile, rc, ri, fc, fi, cp, cp_timestamp):

        # Note: this version of vector_correlate flips the sign of the replica signal rather than the raw signal
        # This was done to validate the method used in CUDARecv.
        # Conceptually, I think they should produce the same result, but this will verify.

        cp_since_prev_bit = (cp-cp_timestamp)%20
        cp_to_next_bit    = 20-cp_since_prev_bit
        idx_next_bit      = np.floor((L_CA*cp_to_next_bit-rc)*(rawfile.fs/fc)).astype(np.int32) + 1

        cp_compl = np.floor((rawfile.S*(fc/rawfile.fs)+rc)/L_CA)

        if (idx_next_bit > 0) and (idx_next_bit < rawfile.S):

            raw = np.array(rawfile.rawsnippet)
            #raw_flip    = np.array(raw_no_flip)
            #raw_flip[idx_next_bit:] = -raw_flip[idx_next_bit:]

            # #### convert to baseband with current fi and ri (velocities)

            doppler_wipeoff_signal = np.exp(-1j*((2.0*PI*fi*rawfile.time_idc) + (2.0*PI*ri)))

            baseband_no_flip = raw*doppler_wipeoff_signal
            #baseband_flip    = raw_flip*doppler_wipeoff_signal

            # #### already in baseband with current fi and ri (velocities)
            # - can do correlation (use shifted replica, will produce code phase error)
            # - can do code wipeoff and fft (will produce carrier freq error)

            # ##### correlation, save code_corr

            r_no_flip = self.chips.take(np.mod(np.floor(rawfile.time_idc*fc + rc),L_CA).astype(np.int32))
            r_flip = np.array(r_no_flip)
            r_flip[idx_next_bit:] = -r_flip[idx_next_bit:]

            rcfft_no_flip = np.conjugate(np.fft.fft(r_no_flip))
            rcfft_flip = np.conjugate(np.fft.fft(r_flip))

            rfft = np.fft.fft(baseband_no_flip)
            corr_no_flip = np.fft.ifft(rcfft_no_flip*rfft)

            corr_flip = np.fft.ifft(rcfft_flip*rfft)

            if np.abs(corr_flip[0])>np.abs(corr_no_flip[0]):
                r = r_flip
                code_corr = np.fft.fftshift(corr_flip)
            else:
                r = r_no_flip
                code_corr = np.fft.fftshift(corr_no_flip)

            # ##### FFT, save carr_fft, carr_fftbins

            baseband = (raw-np.mean(raw)) * r * doppler_wipeoff_signal
            carr_fft = np.fft.fftshift(np.fft.fft(baseband, rawfile.carr_fftpts))

        else:

            raw = np.array(rawfile.rawsnippet)
            doppler_wipeoff_signal = np.exp(-1j*((2.0*PI*fi*rawfile.time_idc) + (2.0*PI*ri)))

            baseband = raw*doppler_wipeoff_signal

            r = self.chips.take(np.mod(np.floor(rawfile.time_idc*fc + rc),L_CA).astype(np.int32))
            rcfft = np.conjugate(np.fft.fft(r))

            rfft = np.fft.fft(baseband)
            corr = np.fft.ifft(rcfft*rfft)

            code_corr    = np.fft.fftshift(corr)

            baseband = (raw-np.mean(raw)) * r * doppler_wipeoff_signal
            carr_fft = np.fft.fftshift(np.fft.fft(baseband, rawfile.carr_fftpts))

        """
        truncSampPrev = np.zeros(rawfile.S, dtype=complex)
        basebandPrev = np.zeros(rawfile.S, dtype=complex)
        truncSampPrev[0:(idx_next_bit)] = raw[0:(idx_next_bit)]
        basebandPrev[0:(idx_next_bit)] = (truncSampPrev[0:(idx_next_bit)] - np.mean(truncSampPrev)) * r[0:(idx_next_bit)] * doppler_wipeoff_signal[0:(idx_next_bit)]
        carr_fftPrev = np.fft.fftshift(np.fft.fft(basebandPrev, rawfile.S))

        truncSampPost = np.zeros(rawfile.S, dtype=complex)
        basebandPost = np.zeros(rawfile.S, dtype=complex)
        truncSampPost[(idx_next_bit):] = raw[(idx_next_bit):]
        basebandPost[(idx_next_bit):] = (truncSampPost[(idx_next_bit):] - np.mean(truncSampPost)) * r[(idx_next_bit):] * doppler_wipeoff_signal[(idx_next_bit):]
        carr_fftPost = np.fft.fftshift(np.fft.fft(basebandPost, rawfile.S))

        carr_fft = carr_fftPrev + carr_fftPost
        """


        while np.max(np.abs(carr_fft)) < np.max(np.abs(code_corr))/2:
            print('WARNING')

            max_code_idx = (np.abs(code_corr).argmax()-rawfile.S/(2.0)).astype(np.int32)
            baseband = (raw-np.mean(raw)) * np.roll(r,max_code_idx) * doppler_wipeoff_signal
            carr_fft = np.fft.fftshift(np.fft.fft(baseband, rawfile.carr_fftpts))


        return code_corr, carr_fft, cp_compl








    def _make_L1_CAcode_chips(self, dtype=np.int32):
        """
        Returns an array of C/A code chips (over values +1 and -1)
        for a given PRN specified in self.prn.

        @type dtype : numpy.dtype
        @param dtype : The data type of the output code chips.
        Default is numpy.int32.
        @rtype : numpy.ndarray
        @return : Array of code chips (over values +1 and -1) for a given PRN
        specified in self.prn.
        """

        prn = self.prn

        assert 1 <= prn <= 210, 'Invalid PRN code:'+str(prn)

        # Initialize the shift register tap points.
        tap1 = np.array([1,0,0,0,0,0,0,1,0,0])
        tap2 = np.array([1,1,1,0,1,0,0,1,1,0])

        # Now initialize the g1 and g2 shift registers.
        g1reg = np.ones(10)#, np.ones(10)
        g2delay,g2reg = self._get_g2_delay_reg ()

        # Initialize g1 and g2 outputs.
        g1 = np.hstack((g1reg, np.zeros(1013)))
        g2 = np.hstack((g2reg, np.zeros(1013)))

        # Generate the g1 and g2 chips based on the feedback polynomials.
        for i in xrange(10, 1023):
            g1[i] = g1reg.dot(tap1) % 2
            g2[i] = g2reg.dot(tap2) % 2
            g1reg = np.append(g1reg[1:10], g1[i])
            g2reg = np.append(g2reg[1:10], g2[i])

        # Shift the g2 code to the right by the delay (in chips).

        g2 = np.roll(g2, g2delay)

        # Form the chips as G1 + G2 modulo 2 taking values of +1 or -1.
        return np.where((g1 + g2) % 2 == 0, -1, 1).astype(dtype)

    def _get_g2_delay_reg (self):
        prn = self.prn
        if prn >= 1 and prn <= 37:
            delays = [\
                5,6,7,8,17,18,139,140,141,251,\
                252,254,255,256,257,258,469,470,471,472,\
                473,474,509,512,513,514,515,516,859,860,\
                861,862,863,950,947,948,950]
            return delays[prn-1],np.ones(10)

        additional_prn = {
            133:  (603,1731),
            135:  (359,1216),
            138:  (386,0450)
        }

        try:
            delay,reg_o = additional_prn[prn]
        except KeyError:
            print 'Delay and G2 register initialization undefined for PRN',prn
            raise

        def oct_to_bin(d):
            s = 0
            p = 1
            while d:
                d, m = divmod(d, 10)
                s = s + m * p
                p = p * 8

            return np.array([np.binary_repr(s,width=10)[k]=='1' for k in range(-10,0)]).astype(int)
        return delay,oct_to_bin(reg_o)

    def _trim_mean(self, arr, percent):
        """
        Returns the trimmed mean using percentage-based trimming.

        @type arr: numpy.ndarray
        @param arr: One-dimensional array of data.
        @type percent: float
        @param percent: 'middle' percent of data that we want to average.
        @rtype: numpy.ndarray
        @return: trimmed mean using percentage-based trimming
        """

        lower = sps.scoreatpercentile(arr, percent/2)
        upper = sps.scoreatpercentile(arr, 100 - percent/2)
        return sps.tmean(arr, limits=(lower, upper), inclusive=(False, False))

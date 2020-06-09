# -*- coding: utf-8 -*-


import filters
import math

class SignalNoiseMeter():
    """
    SignalNoiseMeter instances are objects that can determine the signal to
    noise ratio of a given string of inphase and quadrature prompt correlator
    outputs. (Variance Summing Method)
    """
    
    def __init__(self, N=20, T=0.001): 
        """Initialize this SignalNoiseMeter instance.
        
        @type N : int
        @param N : This is the number of previous samples to average over.
        @type sampletime : float
        @param sampletime : This is the time duration (in seconds) between samples.
        """
            
        self.meanfilter = filters.RunningAverageFilter(N)
        self.varfilter = filters.RunningAverageFilter(N)
        self.averagingtime = N*T
        
    def reset(self, meanpower=0, varpower=0):
        """
        Reset the SignalNoiseMeter.
        
        @type meanpower : float
        @param meanpower : This is the initial state of the running average 
        filter taking power values. The default is zero.
        @type varpower : float
        @param varpower : This is the initial state of the low pass filter 
        taking power variance values. The default is zero.
        """
            
        self.meanfilter.reset(average=meanpower)
        self.varfilter.reset(average=varpower)
        
    def update(self, iP, qP):
        """
        Calculate the current signal to noise ratio using a modification of
        the method presented in SoftGNSS (the Variance Summing Method). The
        difference here is that this function can be called every millisecond
        with new iP and qP inputs because instead of averaging over a fixed
        array of previous iP and qP values, low pass filters are employed to
        provide the mean and variance statistics.
        """
        
        z = iP*iP + qP*qP
        z_mean = self.meanfilter.update(z)
        z_var = self.varfilter.update((z-z_mean)**2)
        sqrtarg = z_mean*z_mean-z_var
        sqrtarg = sqrtarg if sqrtarg > 0 else 0
        carrier_mean = math.sqrt(sqrtarg)
        noise_var = (z_mean - carrier_mean) / 2
        logarg = carrier_mean / (2*self.averagingtime*noise_var)
        logarg = logarg if logarg > 1 else 1
        return 10*math.log10(logarg.__abs__()) # the SNR

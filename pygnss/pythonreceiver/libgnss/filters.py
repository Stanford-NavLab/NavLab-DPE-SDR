# -*- coding: utf-8 -*-


from .constants import *
import collections
import numpy as np

class RunningAverageFilter():
    """
    This is a general purpose filter that outputs the average of the last
    N samples, where N is defined at initialization.
    """
    
    def __init__(self, N, average=0):
        """
        Initialize this filter to output the average of the last N samples.
        
        @type N : int
        @param N : This is the number of previous samples to average over.
        @type average : float
        @param average : This is the initial average value of the initial N samples.
        @rtype : libgnss.RunningAverageFilter
        @return : libgnss.RunningAverageFilter object
        """
        
        self.reset(N=N, average=average)
        
    def reset(self, N=None, average=0):
        """
        Reset the filter history queue by filling it with the desired average value.
        
        @type N : int
        @param N : This is the number of previous samples to average over.  
        When left to the default value of None, N will be set using the current value of N.
        @type average : float
        @param average : This is the desired average value to reset to.
        The default value is zero.
        """
         
        self.N = N if N is not None else self.N
        self.average = average
        self.queue = collections.deque([average]*self.N)
        
    def update(self, xn):
        """
        Returns the average of the last N samples given the most recent sample xn.
        
        @type xn : float
        @param xn : This is the most recent sample.
        @rtype : float
        @return : This is the average of the last N samples, where xn is the Nth sample.
        """
            
        self.average = self.average + (xn - self.queue[0])/self.N
        self.queue[0] = xn
        self.queue.rotate(1)
        return self.average
        
class SingleParameterDigitalFilter():
    """
    This is the parent class of all simple digital filters that only has a single parameter.
    Note: might introducde a DigitalFilter class in the future.
    """
    
    def __init__(self, k):
        """
        Initialize this filter with the parameter k, where k may represent
        anything from a damping ratio to a time interval, depending on its
        usage in child classes.
        
        @type k: float
        @param k: damping ratio or time interval, depending on child classes
        """
        
        self.k = k
        self.reset()
        
    def reset(self, h=0, k=None):
        """
        Reset this filter's history variable h, where h could be anything
        from the most recent output to the accumulated outputs, depending on
        its usage in child classes.
        
        @type h: float
        @param h: Value to reset filter to. The default is 0.
        """
        
        if k is not None:
            self.k = k

        self.h = h
        
class BilinearIntegrator(SingleParameterDigitalFilter):
    """
    Instances of this class are bilinear z-transform integrators.  
    
    For more information, see
     - Misra and Enge, page 478
     - Kaplan and Hegarty, page 181
    """
    
    def update(self, xn):
        """
        Returns the integrated inputs given the most recent input xn.
        
        @type xn: float
        @param xn: input
        @rtype: float
        @return: Integrated inputs, given most recent input xn.
        """   
        # A(n-1) and A(n) are represented here by h0 and self.h, respectively.
        # The time between samples (often denoted T) is represented as self.k.
        h0 = self.h
        self.h = self.h + self.k*xn
        return (self.h + h0) / 2.0
    
class BoxcarIntegrator(SingleParameterDigitalFilter):
    """
    Instances of this class are boxcar z-transform integrators.  
    
    For more information, see
     - Misra and Enge, page 478
     - Kaplan and Hegarty, page 181
    """
        
    def update(self, xn):
        """
        Returns the integrated inputs given the most recent input xn.
        
        @type xn: float
        @param xn: input
        @rtype: float
        @return: Integrated inputs, given most recent input xn.
        """  
        # A(n-1) is represented here by self.h.
        # The time between samples (often denoted T) is represented as self.k.
        self.h = self.h + self.k*xn
        return self.h
        
class LowPassFilter(SingleParameterDigitalFilter):
    """
    The LowPassFilter, unlike the LoopFilter, is a simple general purpose filter
    for smoothing outputs, such as those used in signal to noise measurements and
    lock detection.  
    
    For more information, see
     - Kaplan and Hegarty, page 234, Fig. 5.43
    """
            
    def update(self, xn):
        """
        Returns the filtered inputs given the most recent input xn.
        
        @type xn: float
        @param xn: input
        @rtype: float
        @return: Filtered inputs, given most recent input xn.
        """
        # The most recent output y(n-1) is represented here by self.h.
        self.h = self.k*xn + (1-self.k)*self.h
        return self.h
        
class FIRfilter():
    
    def __init__(self, settings):
        
        self.settings = settings
        
        # remez (FIR) design arguments
        Fs    = self.settings.fs       # sample-rate
        Fpass = 1.5*self.settings.fc   # passband edge
        Fstop = 2.0*self.settings.fc   # stopband edge
        Wp = Fpass/(Fs)    # pass normalized frequency
        Ws = Fstop/(Fs)    # stop normalized frequency
        print(Wp,Ws)
        
        # Create an FIR filter
        # remez takes a list of "bands" and amplitude for each band
        #from scipy.signal import fir_filter_design as ffd
        #self.b = ffd.remez(20, [0, Wp, Ws, .5], [1,0])
        #self.b = np.ones(20)
        #obtained from fft inverse
        self.b = np.array([-0.06331471-0.00048563j, -0.05197719-0.00031893j,
                            0.02641185+0.00012155j,  0.14918167+0.00045769j,
                            0.26136484+0.00040093j,  0.30664062+0.j        ,
                            0.26136484-0.00040093j,  0.14918167-0.00045769j,
                            0.02641185-0.00012155j, -0.05197719+0.00031893j,
                           -0.06331471+0.00048563j])
        self.blen = len(self.b)
        self.prev_array = np.zeros((self.blen-1),dtype=complex)
        
    def update(self,curr_array):
        
        result = np.convolve(self.b,np.concatenate((self.prev_array,curr_array)),'valid')
        self.prev_array = np.array(curr_array[-(self.blen-1):])

        return result
        

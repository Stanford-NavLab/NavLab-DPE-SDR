# -*- coding: utf-8 -*-


import filters

class LockDetector():
    """
    LockDetector instances are objects that can determine whether or not
    a particular tracking channel is locked or not. 
    
    For more information: 
    Kaplan and Hegarty pages 234-235
    
    note to self: Consider implementing false phase lock detection (p. 235) in the 
    LockDetector.update() function.
    """
    
    def __init__(self, N, k, lossthreshold, lockthreshold):
        """
        Initialize the lock detector.
        
        @type N : int
        @param N : This is the number of previous samples to average over.
        @type k : float
        @param k : This is the scale factor (specifically, the divisor of the inphase
        prompt inputs) used for comparing the inphase and quadraphase
        inputs.  A typical value is 1.5.
        @type lossthreshold : int
        @param lossthreshold : Once the loss-of-lock condition has been determined 
        consecutively by the amount given by lossthreshold, loss-of-lock is declared.
        @type lockthreshold : int
        @param lockthreshold : Once the lock condition has been determined consecutively 
        by the amount given by lockthreshold, lock is declared.
        """
        
        self.k = k
        self.lossthreshold = lossthreshold
        self.lockthreshold = lockthreshold
        self.iFilter = filters.LowPassFilter(0.0247)
        self.qFilter = filters.LowPassFilter(0.0247)
        self.reset()
        
    def reset(self, iMagnitude=0, qMagnitude=0, lock=0):
        """
        Reset the lock detector.
        
        @type iMagnitude : float
        @param iMagnitude: The initial state of the inphase running average filter will be set
        to this value (the default is zero).
        @type qMagnitude : float
        @param qMagnitude: The initial state of the quadraphase running average filter will be
        set to this value (the default is zero).
        @type lock : int
        @param lock: The initial state of the lock detector is set to this status (the
        default is 0 for False).
        """
        
        self.losscount = 0
        self.lockcount = 0
        self.iFilter.reset(h=iMagnitude)
        self.qFilter.reset(h=qMagnitude)
        self.lock = lock
        
    def update(self, iP, qP):
        """
        Update the lock detector with the latest prompt correlator outputs,
        and determine the locking status.
        
        @type iP : float
        @param iP : This is the latest inphase prompt correlator output.
        @type qP : float
        @param qP : This is the latest quadraphase prompt correlator output.
        @rtype : tuple
        @return : (int, float), (self.lock, iP-qP)
        This tuple carries locking status (self.lock - either 1 for True or
        0 for False) and the current difference between the filtered
        magnitudes of the (scaled) inphase prompt value and the quadraphase
        prompt value.  The latter can be used as a diagnostic aid when
        reviewing data.
        """
            
        iP = self.iFilter.update(iP.__abs__()) / self.k
        qP = self.qFilter.update(qP.__abs__())
        
        if iP > qP:
            self.losscount = 0
            if self.lockcount > self.lockthreshold:
                self.lock = 1
                
                # Development Notes: consider implementing false phase lock
                # detection logic here.
            else:
                self.lockcount = self.lockcount + 1
        else:
            self.lockcount = 0
            if self.losscount > self.lossthreshold:
                self.lock = 0
            else:
                self.losscount = self.losscount + 1
                
        return self.lock, iP-qP

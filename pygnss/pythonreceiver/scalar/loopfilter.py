# -*- coding: utf-8 -*-


from ..libgnss import filters

class LoopFilter():
    """A channel component for filtering the discriminations."""
    
    def __init__(self, T, channel=None, integrator='BILINEAR', order=2, Bnp=3.0, Bnf=0.0):
        """Constructs loop filter with specified integration period"""
        
        self.T = T
        
        if channel is not None:
            self.channel = channel

        self.integrator = integrator
        self._set_loop_coefficients(order=order, Bnp=Bnp, Bnf=Bnf)                        

    def _set_integration_period(self, T):
        
        self.T = T
        
        if hasattr(self,'intVel') or hasattr(self,'intAcc'):
            self.reset()
            
    def _set_loop_coefficients(self, order = 2, Bnp = 3.0, Bnf = 0.0):
        """Critically damped"""
        
        self.order = order
        
        assert self.order in [2,3], 'Unsupported filter order: '+str(self.order)

        self.Bnp = Bnp
        self.Bnf = Bnf

        if self.order == 2:
            w0p = self.Bnp/0.53
            self.Kvp = w0p**2.0
            self.Kpp = 1.414*w0p
            w0f = self.Bnf/0.25
            self.Kvf = w0f
            if not hasattr(self,'intVel'):
                if self.integrator == 'BOXCAR':
                    self.intVel = filters.BoxcarIntegrator(k=self.T)
                elif self.integrator == 'BILINEAR':
                    self.intVel = filters.BilinearIntegrator(k=self.T)
            if hasattr(self,'intAcc'):
                delattr(self, 'intAcc')
            self.update = self._2nd_order_update
            
        if self.order == 3:
            w0p = self.Bnp/0.7845
            self.Kap = w0p**3.0                        
            self.Kvp = 1.1*w0p**2.0
            self.Kpp = 2.4*w0p
            w0f = self.Bnf/0.53
            self.Kaf = w0f**2.0
            self.Kvf = 1.414*w0f
            if not hasattr(self,'intVel'):
                if self.integrator == 'BOXCAR':
                    self.intVel = filters.BoxcarIntegrator(k=self.T)
                elif self.integrator == 'BILINEAR':
                    self.intVel = filters.BilinearIntegrator(k=self.T)
            if not hasattr(self,'intAcc'):
                if self.integrator == 'BOXCAR':
                    self.intAcc = filters.BoxcarIntegrator(k=self.T)
                elif self.integrator == 'BILINEAR':
                    self.intAcc = filters.BilinearIntegrator(k=self.T)
            self.update = self._3nd_order_update
            
        self.reset()

    def reset(self, intVel = 0, intAcc = 0):
        """
        Reset the accumulation variables in the integrators.
        
        @type  intVel: float
        @param intVel: This is the integrated error in the velocity accumulator. 
                       The default value is zero.
        """
            
        if self.order == 2:
            self.intVel.reset(h=intVel, k=self.T)
        if self.order == 3:
            self.intVel.reset(h=intVel, k=self.T)
            self.intAcc.reset(h=intAcc, k=self.T)
            
    def _2nd_order_update(self, xp = 0, xf = 0):
        """
        Return the filtered discriminator output yn = (xp*Kvp+xf*Kvf)/s + xp*Kpp
        The final integration is a boxcar integration completed outside of the LoopFilter
        
        @type  xp: float
        @param xp: most recent phase discriminator output
        @type  xf: float
        @param xf: most recent frequency discriminator output
        @rtype : float
        @return: filtered output
        """
        
        yn = self.intVel.update(xp*self.Kvp+xf*self.Kvf) + xp*self.Kpp
        
        return yn
        
    def _3nd_order_update(self, xp = 0, xf = 0):
        """
        Return the filtered discriminator output yn = ((xp*Kap+xf*Kaf)/s + xp*Kvp + xf*Kvf)/s + xp*Kpp
        The final integration is a boxcar integration completed outside of the LoopFilter
        
        @type  xp: float
        @param xp: most recent phase discriminator output
        @type  xf: float
        @param xf: most recent frequency discriminator output
        @rtype : float
        @return: filtered output
        """
        
        yn = self.intVel.update( self.intAcc.update(xp*self.Kap+xf*self.Kaf) + xp*self.Kvp + xf*self.Kvf) + xp*self.Kpp
        
        return yn

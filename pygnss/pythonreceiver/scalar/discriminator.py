# -*- coding: utf-8 -*-


from ..libgnss.constants import *
import numpy as np

class Discriminator():
    """A channel component for performing discriminations on correlations."""

    def __init__(self, flavor, channel=None):
        """Constructs Discriminator with specified flavor."""

        if flavor in ['DLL']:
            self.update = self._dll_update
        elif flavor in ['PLL']:
            self.update = self._pll_update
        elif flavor in ['FLL']:
            self.update = self._fll_update
        else:
            print('Warning: unknown discriminator flavor, discriminator update function is not set')

        if channel is not None:
            self.channel = channel

    def _dll_update(self, iE, qE, iL, qL):
        """Outputs DLL discrimination.
        """

        xp = xf = 0.0

        # Use the Normalized Early Minus Late Envelope discriminator.
        E = np.sqrt(iE**2.0 + qE**2.0)
        L = np.sqrt(iL**2.0 + qL**2.0)

        if (E+L) != 0:
            xp = (E-L) / (2.0*(E+L))
            # NOTE: 2.0 because self.offset=0.5 in correlator
        else:
            print("Warning: very low correlation values")

        return xp, xf

    def _pll_update(self, iP, qP):
        """Outputs PLL discrimination.
        """

        xp = xf = 0.0

        if iP != 0:
            xp = np.arctan(qP/iP) / (2.0*PI)
        else:
            print("Warning: very low correlation values")

        return xp, xf

    def _fll_update(self,iP1,qP1,iP0,qP0,N):
        assert (N>1), "Error: N has to be > 1 for FLL, N is "+str(N)
        assert (np.mod(N,2)==0), "Error: N has to be even, N is"+str(N)

        xp = xf = 0.0
        cross = iP0*qP1-iP1*qP0
        dot   = iP0*iP1+qP0*qP1

        if dot > 0.0:
            xf = np.arctan2(cross,dot)/(2.0*np.pi*T_CA* N)
        else:
            xf = np.arctan2(-1.0*cross,-1.0*dot)/(2.0*np.pi*T_CA*N)

        return xp, xf




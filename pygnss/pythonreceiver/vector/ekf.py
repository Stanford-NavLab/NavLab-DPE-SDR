# -*- coding: utf-8 -*-


import numpy as np
from ..libgnss import filters,utils
from ..libgnss.constants import *

class ExtendedKalmanFilter():

    def __init__(self, X_ECEF, Sigma=None, T = 0.020, rxTime_a = None,rxTime0 = None):

        self.X_ECEF = np.matrix(X_ECEF)
        self.Sigma  = Sigma # Note: Sigma = P_k (at various "given k's")
        self.rxTime_a = rxTime_a
        self.rxTime0  = rxTime0

        # initialize container for state correction
        self.dX_ECEF = np.matrix(np.zeros((8,1),dtype=float))

        self.x_list  = []
        self.e_list  = []
        self.e_count = 0
        self.K = np.asmatrix(np.identity(8))
        self._time_update = self._time_update_l5                   # Original
        #self._time_update = self._time_update_m5                    # Debug for CUDARecv
        self._measurement_update = self._measurement_update_l5     # Original
        #self._measurement_update = self._measurement_update_m5      # Debug for CUDARecv

        # Debug for CUDARecv
        self.Sigma = np.asmatrix(np.identity(8))

        # create propagation matrix F for prediction of new state X_{k+1}- from X_{k}+
        tempF      = np.zeros((8,8),dtype=float)
        tempF[0,:] = np.array([1,0,0,0,T,0,0,0])
        tempF[1,:] = np.array([0,1,0,0,0,T,0,0])
        tempF[2,:] = np.array([0,0,1,0,0,0,T,0])
        tempF[3,:] = np.array([0,0,0,1,0,0,0,T])
        tempF[4,:] = np.array([0,0,0,0,1,0,0,0])
        tempF[5,:] = np.array([0,0,0,0,0,1,0,0])
        tempF[6,:] = np.array([0,0,0,0,0,0,1,0])
        tempF[7,:] = np.array([0,0,0,0,0,0,0,1])
        self.F     = np.matrix(tempF)

        # Match CUDARecv by disabling the propagation from velocity to position
        self.F = np.asmatrix(np.identity(8))

        self.H     = np.asmatrix(np.identity(8))

        # initialize container for process noise covariance Q
        self.Q     = np.matrix(np.zeros((8,8),dtype=float))
        self.LPFv  = filters.RunningAverageFilter(20)

    # predict the state X for the next timestamp
    def _predict_X(self):
        self.X_ECEF = self.F*self.X_ECEF
        return self.X_ECEF

    def _update_Q(self):
        # Q is the process noise covariance matrix. Same dimensions as Sigma
        # not too sure of velocity's units
        v = np.linalg.norm(self.X_ECEF[4:7])
        v = self.LPFv.update(v)
        v = 1.0+250.0/(np.min([np.max([v**2.0,50.0]),125.0]))
        Qv = np.matrix(np.zeros((4,4),dtype=float))
        Qv[0,0] = Qv[1,1] = Qv[2,2] = v
        #Qv[3,3] = ((1.0e-7)*(C))**2.0#((2.5e-10)*(C))**2.0
        Qv[3,3] = ((2.5e-10)*(C))**2.0
        Q = np.matrix(np.zeros((8,8),dtype=float))
        Q[4:,4:] = Qv
        self.Q = self.F*Q*(self.F.T)
        #print self.Q
        #self.Q = np.asmatrix(np.identity(8)) # Debug for CUDARecv
        return self.Q

    def _predict_Sigma(self):
        # Compute P_k|k1 = F_k*P_k1|k1*F^T_k + Q_k
        self.Sigma = self.F*self.Sigma*(self.F.T) + self.Q
        return self.Sigma

    def _update_K(self, H, W):
        # Compute S^-1_k = inv(H_k*P_k|k1*H_k^T + R_k)
        S_inv = ((H*self.Sigma*(H.T)+W).I)
        # Compute K_k = P_k|k1*H^T_k*S^-1_k
        self.K = self.Sigma*(H.T)*S_inv
        return self.K

    def _update_dX(self, e):
        # Compute K_k*y_k
        self.dX_ECEF = self.K*e
        return self.dX_ECEF

    def _correct_X(self):
        # Compute x_k|k = x_k|k1 + K_k*y_k
        self.X_ECEF = self.X_ECEF + self.dX_ECEF
        return self.X_ECEF

    def _update_dX_ECI(self, e):
        # Compute K_k*y_k
        self.dX_ECI = self.K*e
        #dX_ECI[3,0] = dX_ECI[3,0]/c
        #dX_ECI[7,0] = dX_ECI[7,0]/c
        #self.dX_ECI = dX_ECI
        return self.dX_ECI

    def _correct_X_ECI(self):
        # Compute x_k|k = x_k|k1 + K_k*y_k
        self.rxTime_a = self.rxTime0 - (self.X_ECEF[3,0] + self.dX_ECI[3,0])/C
        X_ECI = utils.ECEF_to_ECI(self.X_ECEF, t_gps = self.rxTime_a, t_c = self.rxTime_a)
        X_ECI = X_ECI + self.dX_ECI
        self.X_ECEF = utils.ECI_to_ECEF(X_ECI, t_gps = self.rxTime_a, t_c = self.rxTime_a)
        return self.X_ECEF, X_ECI

    def _correct_Sigma(self, H):
        # Compute P_k|k = (I-K*H)*P_k|k1
        self.Sigma = (np.asmatrix(np.identity(8)) - self.K*H)*self.Sigma
        return self.Sigma

    def _update_W(self, e):
        # Note: W here is R_k in Wikipedia Kalman equations form
        # use a queue
        self.e_list.append(e)
        W = np.matrix(np.cov(np.transpose(np.asarray(self.e_list)[:,:,0])))
        self.e_list.pop(0)
        return W/np.sqrt(2)

    def _init_Sigma(self):
        self.Sigma = np.matrix(np.cov(np.transpose(np.asarray(self.x_list)[:,:,0])))
        return self.Sigma

    def _measurement_update(self):
        pass

    def _time_update(self):
        pass

    def _measurement_update_l5(self, e):

        dX_ECEF = self._update_dX(e)
        X_ECEF  = self._correct_X()

        '''
        self.x_list.append(X_ECEF)
        self.e_list.append(e)
        self.e_count = self.e_count+1
        if self.e_count > 5:
            self._time_update = self._time_update_m5
            self._measurement_update = self._measurement_update_m5
            self._init_Sigma()
        return
        '''
        return

    def _time_update_l5(self):

        X_ECEF = self._predict_X()

        return

    def _measurement_update_m5(self, e):

        #W       = self._update_W(e) # Debug for CUDARecv; uncomment for normal PyGNSS
        W       = np.asmatrix(np.identity(8))
        K       = self._update_K(self.H,W)
        dX_ECEF = self._update_dX(e)
        X_ECEF  = self._correct_X()
        Sigma   = self._correct_Sigma(self.H)

        return e, dX_ECEF, X_ECEF, Sigma

    def _time_update_m5(self):

        # PREDICT: new X, new Sigma
        X_ECEF = self._predict_X()
        Q      = self._update_Q()
        Sigma  = self._predict_Sigma()

        return X_ECEF, self.Q, Sigma

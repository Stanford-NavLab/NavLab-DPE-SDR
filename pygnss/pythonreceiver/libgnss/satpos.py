# -*- coding: utf-8 -*-

from constants import *

import utils
import numpy as np

def locate_satellite(ephemerides, ctime, clkb, clkd):
    """
    Return the state (position, clock bias, velocity, and clock drift) 
    of the requested satellite at the requested corrected GPS transmit time.
    The satellite state is calculated using the ephemerides. 
    
    For more information, see
    Kaplan and Hergarty page 42 for position calculation from ephemerides.
    Computing satellite velocity using the broadcast ephemeris by Benjamin W. Remondi
    for velocity calculation from ephemerides        
     - http://fenrir.naruoka.org/download/autopilot/note/080205_gps/gps_velocity.pdf
     - http://www.ngs.noaa.gov/gps-toolbox/bc_velo/bc_velo.c        
    
    @type  ephemerides: libgnss.ephemerides.Ephemerides
    @param ephemerides: Satellite orbital parameters used to calculate the state 
    (position, clock bias, velocity, clock drift) of the specific satellite. 
    @type  ctime: float
    @param ctime: Corrected GPS transmit time. This is the time in GPS 
    seconds of week for which we want to compute the satellite state. 
    @type  clkb: float
    @param clkb: Clock bias of the satellite, used to populate state matrix.
    @type  clkd: float
    @param clkd: Clock drift of the satellite, used to populate state matrix.
    @rtype : numpy.matrix
    @return: Satellite state (position, clockbias, velocity, clock drift) 
    of the requested satellite at the requested corrected GPS transmit time.
    """
        
    # Establish a shorthand to make the code less verbose.
    eph = ephemerides
    # Initialize output matrix
    stateMat = np.matrix(np.zeros((8,1)))
    
    # Position calculation matches page 42 of Kaplan and Hegarty 
    
    # Semi-major axis
    a = eph.sqrt_A**2.0
    # Corrected mean motion
    n = np.sqrt(MU / (a**3.0)) + eph.delta_n
    
    # Compute satellite position and velocity
    tk = correct_week_crossover(ctime - eph.t_oe)
    # Mean anomaly
    E = M = np.mod(eph.M_0 + n*tk, 2.0*PI)
    # Eccentric anomaly
    for count in range(10):
        f = M - E + eph.e * np.sin(E)
        dfdE = - 1.0 + eph.e * np.cos(E)
        dE = - f / dfdE
        E = np.mod(E + dE, 2.0*PI)
        if (np.abs(dE) < 1.0e-12):
            # Convergence is achieved, break the loop prematurely.
            break
    if (np.abs(dE) > 1.0e-12):
        print('Warning: Iterative compution of the eccentric anomaly did not converge.')
    sinE = np.sin(E)
    cosE = np.cos(E)
    
    # True anomaly
    v = np.arctan2(np.sqrt(1.0-eph.e**2.0)*sinE/(1.0-eph.e*cosE), (cosE-eph.e)/(1.0-eph.e*cosE))
    # Argument of latitude
    u = np.mod(v + eph.omega, 2.0*PI)
    
    # second harmonic perturbations
    cos2u = np.cos(2.0*u)
    sin2u = np.sin(2.0*u)
    
    # argument of latitude correction  
    d_u = eph.C_uc * cos2u + eph.C_us * sin2u
    # radius correction  
    d_r = eph.C_rc * cos2u + eph.C_rs * sin2u
    # orbital inclination correction
    d_i = eph.C_ic * cos2u + eph.C_is * sin2u
    
    # corrected argument of latitude
    u = u + d_u
    # corrected radius
    r = a * (1.0 - eph.e * cosE) + d_r
    # corrected inclination
    i = eph.i_0 + eph.IDOT * tk + d_i        
    # corrected longitude of node
    omegak = np.mod(eph.OMEGA_0 + (eph.OMEGADOT - OEDot)*tk - OEDot*(eph.t_oe),2.0*PI)

    # positions in orbital plane
    x_op = r * np.cos(u)
    y_op = r * np.sin(u)

    cos_omegak = np.cos(omegak)
    sin_omegak = np.sin(omegak)
    cosi = np.cos(i)
    sini = np.sin(i)
    
    stateMat[0,0] = x_op * cos_omegak - y_op * sin_omegak * cosi
    stateMat[1,0] = x_op * sin_omegak + y_op * cos_omegak * cosi
    stateMat[2,0] = y_op * sini
    stateMat[3,0] = clkb
    
    # Velocity calculation matches bc_velo.c given by NGS NOAA

    # second harmonic perturbations
    cos2u = np.cos(2.0*u)
    sin2u = np.sin(2.0*u)

    edot  = n / (1.0 - eph.e*cosE)
    vdot  = sinE*edot*(1.0 + eph.e*np.cos(v)) / (np.sin(v)*(1.0-eph.e*cosE))  
    udot  = vdot + 2.0*(eph.C_us*cos2u - eph.C_uc*sin2u)*vdot
    rdot  = a*eph.e*sinE*edot + 2.0*(eph.C_rs*cos2u - eph.C_rc*sin2u)*vdot
    idotdot = eph.IDOT + (eph.C_is*cos2u - eph.C_ic*sin2u)*2.0*vdot

    vx_op = rdot*np.cos(u) - y_op*udot
    vy_op = rdot*np.sin(u) + x_op*udot
    omegadot = eph.OMEGADOT - OEDot

    tmpa = vx_op - y_op*cosi*omegadot  
    tmpb = x_op*omegadot + vy_op*cosi - y_op*sini*idotdot

    stateMat[4,0] = tmpa * cos_omegak - tmpb * sin_omegak;  
    stateMat[5,0] = tmpa * sin_omegak + tmpb * cos_omegak;  
    stateMat[6,0] = vy_op*sini + y_op*cosi*idotdot;                
    
    stateMat[7,0] = clkd
    
    return stateMat

def satellite_clock_correction(ephemerides, transmitTime):
    """
    Helper function that returns the clock correction (clock bias, clock drift) 
    of the requested satellite at the requested transmit time. The satellite 
    clock correction is calculated using the ephemerides decoded in subframe 1. 
    
    note: values are returned in seconds and seconds/seconds
    
    ICD-GPS-200C, p. 88-101 
    line 00072 http://gnsstk.sourceforge.net/gps_8c-source.html
    
    @type  ephemerides: libgnss.ephemerides.Ephemerides
    @param ephemerides: Satellite orbital parameters used to calculate the state 
    (in this case, clock correction using subframe 1) of the specific satellite. 
    @type  transmitTime: float
    @param transmitTime: This is the satellite transmit time in GPS seconds 
    of week for which we want to compute the satellite clock correction.
    @rtype: tuple
    @return: (clockbias, clockdrift) in seconds and seconds/seconds
    """
    
    # Establish a shorthand to make the code less verbose.
    eph = ephemerides
    time = transmitTime
    
    # Semi-major axis
    a = eph.sqrt_A**2.0
    # Corrected mean motion
    n = np.sqrt(MU / (a**3.0)) + eph.delta_n
    
    # Compute satellite clock corrections
    tc = correct_week_crossover(time - eph.t_oc) # without corrections
    clkb = eph.a_f2*tc*tc + eph.a_f1*tc + eph.a_f0 - eph.T_GD # without relativistic correction
    tk = correct_week_crossover(time - clkb - eph.t_oe) # without relativistic correction
    # Mean anomaly
    E = M = np.mod(eph.M_0 + n*tk, 2.0*PI)
    # Eccentric anomaly
    for count in range(10):
        f = M - E + eph.e * np.sin(E)
        dfdE = - 1.0 + eph.e * np.cos(E)
        dE = - f / dfdE
        E = np.mod(E + dE, 2.0*PI)
        if (np.abs(dE) < 1.0e-12):
            # Convergence is achieved, break the loop prematurely.
            break
    if (np.abs(dE) > 1.0e-12):
        print('Warning: Iterative compution of the eccentric anomaly did not converge.')
    dtr = F*eph.e*eph.sqrt_A*np.sin(E)
    #tc = correct_week_crossover(time - (clkb + dtr) - eph.t_oc)
    tc = time - (clkb + dtr) - eph.t_oc
    clkb = eph.a_f2*tc*tc + eph.a_f1*tc + eph.a_f0 + dtr - eph.T_GD
    clkd = eph.a_f1 + 2.0*eph.a_f2*tc

    return clkb, clkd
    
def correct_week_crossover(times):
    """
    Helper function that corrects the given times for the half week crossover.
    
    @type  times : numpy.ndarray or int
    @param times : Array of times (in seconds) to correct.
    @rtype : numpy.ndarray or int
    @return: Array of corrected times (in seconds).
    """
        
    t = np.where(times > 302400.0, times - 604800.0, times)
    return np.where(t < -302400.0, t + 604800.0, t)
'''   
iono = {}
iono['a0'] = .1676e-07
iono['a1'] = .2235e-07
iono['a2'] = -.1192e-06
iono['a3'] = -.1192e-06
iono['b0'] = .1106e+06
iono['b1'] = .9830e+05
iono['b2'] = -.1311e+06
iono['b3'] = -.1966e+06

ephemerides, transmitTime, X_ECI, sat_ECI, iono

def ionospheric_correction(ephemerides, transmitTime, X_ECI, sat_ECI, iono):
    #calculate ionospheric delay
    #based on ionocorr() from http://home-2.worldonline.nl/~samsvl/stdalone.pas
    
    from math import radians, cos, sin
    
    llh = utils.geodetic()
    
    # convert to semi-circles
    Latu = radians(llh.lat) / pi
    Lonu = radians(llh.lon) / pi
    Az = radians(satinfo.azimuth[svid]) / pi
    El = radians(satinfo.elevation[svid]) / pi

    a0 = satinfo.ionospheric[svid].a0
    a1 = satinfo.ionospheric[svid].a1
    a2 = satinfo.ionospheric[svid].a2
    a3 = satinfo.ionospheric[svid].a3
    b0 = satinfo.ionospheric[svid].b0
    b1 = satinfo.ionospheric[svid].b1
    b2 = satinfo.ionospheric[svid].b2
    b3 = satinfo.ionospheric[svid].b3
    
    # main calculation
    phi = 0.0137 / (El + 0.11) - 0.022
    Lati = Latu + phi * cos(Az * pi)
    if Lati > 0.416:
        Lati = 0.416
    elif Lati < -0.416:
        Lati = -0.416
        
    Loni = Lonu + phi * sin(Az * pi) / cos(Lati * pi)
    Latm = Lati + 0.064 * cos((Loni - 1.617) * pi)
    
    T = 4.32E+4 * Loni + transmitTime
    
    while T >= 86400:
        T = T - 86400
    while T < 0:
        T = T + 86400
        
    F = 1.0 + 16.0 * (0.53 - El) * (0.53 - El) * (0.53 - El)
    
    per = b0 + b1 * Latm + b2 * Latm * Latm + b3 * Latm * Latm * Latm
    
    if per < 72000.0:
        per = 72000.0
    x = 2 * pi * (T - 50400.0) / per
    amp = a0 + a1 * Latm + a2 * Latm * Latm + a3 * Latm * Latm * Latm
    if amp < 0.0:
        amp = 0.0
    if abs(x) >= 1.57:
        dTiono = F * 5.0E-9
    else:
        dTiono = F * (5.0E-9 + amp * (1.0 - x * x / 2.0 + x * x * x * x /24.0))
    return dTiono * c #m

def tropospheric_correction_standard(satinfo, svid):
    #tropospheric correction using standard atmosphere values
    from math import sin, sqrt, radians

    El = radians(satinfo.elevation[svid])
    
    dRtrop = 2.312 / sin(sqrt(El * El + 1.904E-3)) + 0.084 / sin(sqrt(El * El + 0.6854E-3))
    return dRtrop #m
'''

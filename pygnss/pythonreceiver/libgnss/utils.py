# -*- coding: utf-8 -*-


from .constants import *
import numpy as np

_IE24  = {'a': 6378388.0, 'invf': 297.0}
_IE67  = {'a': 6378160.0, 'invf': 298.247}
_WGS72 = {'a': 6378135.0, 'invf': 298.26}
_GRS80 = {'a': 6378137.0, 'invf': 298.257222101}
_WGS84 = {'a': 6378137.0, 'invf': 298.257223563}

def ECEF_to_LLA(posvel_ECEF, ellipsoid=_WGS84, normalize=False, in_degrees=True):
    """
    Returns lla position in a record array with keys 'lat', 'lon' and 'alt'.
    lla is calculated using the closed-form solution.

    For more information, see
     - Datum Transformations of GPS Positions
       https://microem.ru/files/2012/08/GPS.G1-X-00006.pdf
       provides both the iterative and closed-form solution.

     >>> #ECE Building and Mount Everest
     >>> array([(40.11497089608554, -88.22793631642435, 203.9925799164921),
                (27.98805865809616, 86.92527453636706, 8847.923165871762)],
         dtype=[('lat', '<f8'), ('lon', '<f8'), ('alt', '<f8')])

    @type  posvel_ECEF: np.ndarray
    @param posvel_ECEF: ECEF position in array of shape (3,N)
    @type  ellipsoid: dict
    @param ellipsoid: Reference ellipsoid (eg. _WGS84 = {'a': 6378137.0, 'invf': 298.257223563}).
    Class headers contains some pre-defined ellipsoids: _IE24, _IE67, _WGS72, _GRS80, _WGS84
    @type  normalize: bool
    @param normalize: Its default value is False; True will cause longitudes returned
    by this function to be in the range of [0,360) instead of (-180,180].
    Setting to False will cause longitudes to be returned in the range of (-180,180].
    @rtype : numpy.ndarray
    @return: Position in a record array with keys 'lat', 'lon' and 'alt'.
    """

    from numpy import sqrt, cos, sin, pi, arctan2 as atan2
    #python subtlety only length 1 arrays can be converted to Python scalars
    #if we use the functions from math, we have to use map

    a    = ellipsoid['a']
    invf = ellipsoid['invf']
    f = 1.0/invf
    b = a*(1.0-f)
    e  = sqrt((a**2.0-b**2.0)/a**2.0)
    ep = sqrt((a**2.0-b**2.0)/b**2.0)

    xyz = np.asarray(posvel_ECEF)
    x = xyz[0,:]
    y = xyz[1,:]
    z = xyz[2,:]

    # Create the record array.
    cols = np.shape(xyz)[1]
    lla  = np.zeros(cols, dtype = { 'names' : ['lat', 'lon', 'alt'],
                                  'formats' : ['<f8', '<f8', '<f8']})

    lon = atan2(y, x)
    p = sqrt(x**2.0+y**2.0)
    theta = atan2(z*a,p*b)
    lat = atan2((z+(ep**2.0)*(b)*(sin(theta)**3.0)),(p-(e**2.0)*(a)*(cos(theta)**3.0)))
    N = a/sqrt(1.0-((e**2.0)*(sin(lat)**2.0)))
    alt = p/cos(lat)-N

    if normalize:
        lon = np.where(lon < 0.0, lon + 2.0*pi, lon)

    lla['alt'] = alt

    if in_degrees:
        lla['lat'] = lat*180.0/pi
        lla['lon'] = lon*180.0/pi
    else:
        lla['lat'] = lat
        lla['lon'] = lon

    return lla

def LLA_to_ECEF(lat,lon,alt, ellipsoid=_WGS84):
    '''
    For more information, see
     - Datum Transformations of GPS Positions
       https://microem.ru/files/2012/08/GPS.G1-X-00006.pdf
       provides both the iterative and closed-form solution.

    '''
    from numpy import sqrt, cos, sin, pi, arctan2 as atan2
    #python subtlety only length 1 arrays can be converted to Python scalars
    #if we use the functions from math, we have to use map

    a    = ellipsoid['a']
    invf = ellipsoid['invf']
    f = 1.0/invf
    b = a*(1.0-f)
    e  = sqrt((a**2.0-b**2.0)/a**2.0)

    lat = np.deg2rad(lat)
    lon = np.deg2rad(lon)
    alt = np.array(alt)

    N = a/sqrt(1-e**2*sin(lat)**2)
    X = (N+alt)*cos(lat)*cos(lon)
    Y = (N+alt)*cos(lat)*sin(lon)
    Z = (((b**2)/(a**2))*N+alt)*sin(lat)

    posvel_ECEF = np.matrix(np.zeros([3,len(lat)]))
    posvel_ECEF[0,:] = X
    posvel_ECEF[1,:] = Y
    posvel_ECEF[2,:] = Z

    return posvel_ECEF

def ECI_to_ECEF(posvel_ECI, t_gps=None, t_c=None):
    """
    Returns rotated ECEF position and velocity coordinates
    due to the Earth's rotation during the given time duration from t_c to t_gps.

    For more information, see
     - Convert ECI to ECEF Coordinates by Darin Koblick
       http://www.mathworks.com/matlabcentral/fileexchange/28233-convert-eci-to-ecef-coordinates/content/ECI2ECEF.zip
     - Constell user manual for velocity rotation.
       http://www.constell.org/Downloads/Manual-Version70.pdf
       (not sure if it is correct,
       there is a GPS World article that verified Constell)
     - NIST Technical Note 1385
       Global Position System Receivers and Relativity

    @type  posvel_ECI: np.matrix
    @param posvel_ECI: (8,N) np.matrix
    @type  t_gps: float
    @param t_gps: This is the time coordinate (s) of the rotating earth-fixed frame (ECEF).
    @type  t_c: float
    @param t_c: This is the time coordinate (s) of the non-rotating inertial reference frame (ECI).
    @rtype : np.matrix
    @return: (8,N) np.matrix with rotated position and velocity.
    """

    xyz    = posvel_ECI[0:3,:]
    xyzdot = posvel_ECI[4:7,:]

    # Compute the rotation angle.
    otau = OEDot*(t_gps-t_c)

    # Establish the rotation matrix.
    cotau  = np.cos(otau)
    sotau  = np.sin(otau)
    rot    = np.matrix([[cotau, sotau, 0], [-sotau, cotau, 0], [0, 0, 1]])
    rotdot = np.matrix([[0, - OEDot, 0], [OEDot, 0, 0], [0, 0, 0]])

    # Perform the rotation.
    rotxyz    = rot*xyz
    rotxyzdot = rot*(xyzdot - rotdot*xyz)

    #the following equations should give the same result
    #cotau_oedot = np.cos(otau)*OEDot
    #sotau_oedot = np.sin(otau)*OEDot
    #rotdot = np.matrix([[-sotau_oedot, cotau_oedot, 0], [-cotau_oedot, -sotau_oedot, 0], [0, 0, 0]])
    #rotxyzdot = rot*xyzdot + rotdot*xyz

    posvel_ECEF = np.matrix(np.zeros(np.shape(posvel_ECI)))
    posvel_ECEF[0:3,:] = rotxyz
    posvel_ECEF[3,:]   = posvel_ECI[3,:]
    posvel_ECEF[4:7,:] = rotxyzdot
    posvel_ECEF[7,:]   = posvel_ECI[7,:]

    return posvel_ECEF


def ECEF_to_ECI(posvel_ECEF, t_gps=None, t_c=None):
    """
    Returns rotated ECI position and velocity coordinates
    due to the Earth's rotation during the given time duration from t_c to t_gps.

    For more information, see
     - IS-GPS-200H page 106 for position rotation.
     - Convert ECI to ECEF Coordinates by Darin Koblick
       http://www.mathworks.com/matlabcentral/fileexchange/28233-convert-eci-to-ecef-coordinates/content/ECI2ECEF.zip
     - Constell user manual for velocity rotation.
       http://www.constell.org/Downloads/Manual-Version70.pdf
       (not sure if it is correct,
       there is a GPS World article that verified Constell)
     - NIST Technical Note 7385
       Global Position System Receivers and Relativity
     - "Fundamentals of Inertial Navigation, Satellite-based
       Positioning and their Integration" by A. Noureldin et al,
       DOI: 10.1007/978-3-642-30466-8_2 for coordinate transformations

    @type  posvel_ECEF: np.matrix
    @param posvel_ECEF: (8,N) np.matrix
    @type  t_gps: float
    @param t_gps: This is the time coordinate (s) of the rotating earth-fixed frame (ECEF).
    @type  t_c: float
    @param t_c: This is the time coordinate (s) of the non-rotating inertial reference frame (ECI).
    @rtype : np.matrix
    @return: (8,N) np.matrix with rotated position and velocity.
    """

    xyz    = posvel_ECEF[0:3,:]
    xyzdot = posvel_ECEF[4:7,:]

    # Compute the rotation angle.
    otau = OEDot*(t_gps-t_c)

    # Establish the rotation matrix.
    cotau = np.cos(otau)
    sotau = np.sin(otau)
    rot = np.matrix([[cotau, -sotau, 0], [sotau, cotau, 0], [0, 0, 1]])
    rotdot = np.matrix([[0, -OEDot, 0], [OEDot, 0, 0], [0, 0, 0]])
    # Perform the rotation.
    rotxyz = rot*xyz
    rotxyzdot = rot*xyzdot + rotdot*rotxyz
    #following equations should give the same result
    #cotau_oedot = np.cos(otau)*OEDot
    #sotau_oedot = np.sin(otau)*OEDot
    #rotdot = np.matrix([[-sotau_oedot, -cotau_oedot, 0], [cotau_oedot, -sotau_oedot, 0], [0, 0, 0]])
    #rotxyzdot = rot*xyzdot + rotdot*xyz

    posvel_ECI = np.matrix(np.zeros(np.shape(posvel_ECEF)))
    posvel_ECI[0:3,:] = rotxyz
    posvel_ECI[3,:]   = posvel_ECEF[3,:]
    posvel_ECI[4:7,:] = rotxyzdot
    posvel_ECI[7,:]   = posvel_ECEF[7,:]

    return posvel_ECI

def ECEF_to_ENU(refState=None, curState=None, diffState=None):
    """
    Returns 3xN matrix in enu order
    - Toolbox for attitude determination
       Zhen Dai, ZESS, University of Siegen, Germany

    @type  curState: 3x1 or 8x1 matrix
    @param curState: current position in ECEF
    @type  refState: 3x1 or 8x1 matrix
    @param refState: reference position in ECEF
    @rtype : tuple
    @return: (numpy.matrix (3,N) (ENU), numpy.matrix (3,3) R_ECEF2ENU)
    """

    if refState is None:
        refState = np.matrix([[151055.3983],[-4882530.31559],[4087649.46970]])

    lla = ECEF_to_LLA(refState,in_degrees=False)
    lat = lla['lat'][0]
    lon = lla['lon'][0]

    slon = np.sin(lon)
    clon = np.cos(lon)
    slat = np.sin(lat)
    clat = np.cos(lat)

    R_ECEF2ENU = np.mat([[ -slon, clon, 0.0],
                         [ -slat*clon, -slat*slon, clat],
                         [  clat*clon,  clat*slon, slat]])

    xyz0 = refState[0:3,:]

    if (curState is not None) and (diffState is None):
        xyz  = curState[0:3,:]
        dxyz = xyz - np.tile(xyz0, (1,np.shape(curState)[1]))
    elif (curState is None) and (diffState is not None):
        dxyz = diffState[0:3,:]
    else:
        print("Unknown error, shape of states:")
        print('refState: '+str(np.shape(refState)))
        print('curState: '+str(np.shape(curState)))
        print('diffState: '+str(np.shape(diffState)))

    matENU = R_ECEF2ENU*dxyz

    return matENU, R_ECEF2ENU

def ENU_to_ECEF(refState=None, diffState=None, R_ECEF2ENU= None):
    """
    """
    if R_ECEF2ENU is None:
        lla = ECEF_to_LLA(refState, in_degrees=False)
        lat = lla['lat'][0]
        lon = lla['lon'][0]

        slon = np.sin(lon)
        clon = np.cos(lon)
        slat = np.sin(lat)
        clat = np.cos(lat)

        R_ECEF2ENU = np.mat([[ -slon, clon, 0.0],
                             [ -slat*clon, -slat*slon, clat],
                             [  clat*clon,  clat*slon, slat]])

    R_ENU2ECEF = R_ECEF2ENU.T

    xyz0 = refState[0:3,:]
    dxyz = diffState[0:3,:]

    matECEF = R_ENU2ECEF*dxyz+np.tile(xyz0, (1,np.shape(dxyz)[1]))

    return matECEF

def ENU_to_elaz(ENU):

    ENU = np.asarray(ENU)
    east  = ENU[0,:]
    north = ENU[1,:]
    up    = ENU[2,:]

    # Create the record array.
    cols = np.shape(ENU)[1]
    elazd = np.zeros(cols, dtype = { 'names' : ['ele', 'azi', 'dist'],
                                   'formats' : ['<f8', '<f8', '<f8']})

    horz_dist = np.sqrt(east**2 + north**2)
    elazd['ele']  = np.arctan2(up, horz_dist)
    elazd['azi']  = np.arctan2(east, north)
    elazd['dist'] = np.sqrt(east**2 + north**2 + up**2)

    return elazd

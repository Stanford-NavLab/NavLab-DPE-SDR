import datetime
import pytz
from pythonreceiver import receiver
from pythonreceiver.scalar import channel
from pythonreceiver.libgnss import ephemeris,utils
import numpy as np


FT_PER_M = 1 #100, in meters, changed by Shubhendra /(2.54 * 12)
# time and index
# position in ECEF, LLA, XYZ (per selected ref)

def header(f):#:, euler = False):
    f.write('{0:>7}, {1:>8}, {2:>13},{3:>14}X,{3:>14}Y,{3:>14}Z,{4:>8}X,{4:>8}Y,{4:>8}Z,'.format(
            'Count#','Date','Time','WGS84_p','WGS84_v'
    ))

    f.write('{:>12},{:>12},{:>12},'.format('Lat','Lon','Alt'))

    #f.write('{0:>14}X,{0:>14}Y,{0:>14}Z,{1:>8}X,{1:>8}Y,{1:>8}Z,{2:>14}X,{2:>14}Y,{2:>14}Z,'.format(
    #           'ENU_p','ENU_v','Ref_WGS84_p'
    #))
    f.write('{:>9},{:>9},{:>9}'.format('Heading','Pitch','Roll'))
    #f.write(', SV')
    f.write('\n')

def printer (mc,weekno,rxTime_a,posvel_ECEF,f,hdg = None, pitch = None, roll = None):
    # time, index, position in ECEF, LLA, XYZ (ft, deg)
    # velocity in ECEF, LLA,  XYZ
    # attitude (heading, pitch, roll)
    # position of rxf.write('\n')0
    # velocity of rx0
    # ...

    #receiver = rx

    #rxTime_a, rxTime, posvel_ECEF, posvel_ECI, sats_ECI = naveng.calculate_nav_soln(receiver,mc=mc)
    '''
    prn = list(rx.channels.keys())
    for k in rx.channels.keys():
        try:
            rx.channels[k].ephemerides
        except AttributeError:
           print 'Channel',k,'has no ephemeris'
        prn.remove(k)
    '''
    gpst = datetime.datetime(1980,1,6,0,0,0,0,pytz.utc) + datetime.timedelta(days = weekno * 7,seconds = rxTime_a-18)

    f.write('{:7d}, '.format(mc)) # count
    f.write(gpst.strftime('%Y%m%d, %H%M%S.%f,')) # time (8+1+13)

    ecef_list = (np.array(posvel_ECEF)*FT_PER_M).flatten().tolist()
    ecef_list = np.reshape(ecef_list,(8,)).tolist()
    #print ecef_list
    f.write(('%+15.3f,'*3)%tuple(ecef_list[ :3]))
    f.write(('%+9.3f,' *3)%tuple(ecef_list[4:7]))

    lla_list = utils.ECEF_to_LLA(posvel_ECEF)[0]
    f.write(('%+12.6f,%+12.6f,%+12.3f,')%(lla_list['lat'],lla_list['lon'],lla_list['alt']*FT_PER_M))

    #enu_p,R = utils.ECEF_to_ENU(refState = ref_ECEF,curState = posvel_ECEF)
    #enu_v   = R * posvel_ECEF [4:7]

    #f.write(('%+15.3f,'*3)%tuple((enu_p*FT_PER_M).flatten().tolist()[0]))
    #f.write(('%+9.3f,' *3)%tuple((enu_v*FT_PER_M).flatten().tolist()[0]))
    #f.write(('%+15.3f,'*3)%tuple((ref_ECEF*FT_PER_M).flatten().tolist()[0][:3]))
    f.write('%9.3f,'% (hdg if hdg is not None else 0.0))
    f.write('%+9.3f,'% (pitch if pitch is not None else 0.0))
    f.write('%+9.3f'% (roll if roll is not None else 0.0))
    #f.write(', %2d'%sv)
    f.write('\n')

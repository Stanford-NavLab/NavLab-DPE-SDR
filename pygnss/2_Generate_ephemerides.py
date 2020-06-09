
execfile('setting.py')

from pythonreceiver.libgnss.constants import *
from pythonreceiver.libgnss import rawfile,utils,satpos
from pythonreceiver.scalar import channel, correlator, discriminator, loopfilter,naveng
from pythonreceiver import receiver

import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio

#import printer,os
#from pythonreceiver.libgnss import pygmaps
from pythonreceiver.vector.ekf import ExtendedKalmanFilter

usrp = []

for ip in ip_list:
    rfile = rawfile.RawFile(
        metafile= None,
        abspath = datpath + prefix[:15]+'_usrp'+str(ip)+'_%dkHz.dat'%(fs/1e3),
        fs = fs, fi = 0.0e6, ds = 1.0,
        datatype = np.dtype([('i', np.short), ('q', np.short)]),
        notes = 'Data set '+ prefix[:-1]
    )

    usrp += [receiver.Receiver(rfile)]

    rx = usrp[-1]
    rx.load_measurement_logs(dirname = prepath, subdir= 'end-of-%d_usrp'%proc_time+ str(ip))
    try:
        os.makedirs(prepath + 'eph%d'%ip)
    except:
        pass

    for prn in rx.channels:
        try:
            rx.parse_ephemerides(prn_list = [prn],m_start=40)
            rx.channels[prn].ephemerides.save_ephemerides(prepath+'eph%d/channel%d.mat'%(ip,prn))
        except:
            pass



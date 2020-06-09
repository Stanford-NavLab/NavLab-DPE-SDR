
execfile('setting.py')

### Main code starts

from pythonreceiver.libgnss.constants import *
from pythonreceiver.libgnss import rawfile,utils,satpos,ephemeris
from pythonreceiver.scalar import channel, correlator, discriminator, loopfilter,naveng
from pythonreceiver import receiver
import printer
import threading

import numpy as np
import scipy.io as sio
import time,os

scalar_usrp = [
    receiver.Receiver(\
        rawfile.RawFile(
            metafile = None,
            abspath  = datpath + refname + prefix[:15] + '_usrp'+str(ip)+'_%dkHz.dat'%int(fs/1e3),
            fs = 2.5e6, fi = 0.0e6, ds = 1.0,
            datatype = np.dtype([('i', np.short), ('q', np.short)]),
            notes = 'Data set '+ refname + prefix[:15]
        )
        , mcount_max = run_time * 1000 + 10000
    ) for ip in ip_list
]

print 'Start scalar tracking @',init_time,'for',run_time

class scalar_thread (threading.Thread):
    def __init__(self,rx,ip,lock):
        threading.Thread.__init__(self)
        self.rx = rx
        self.ip = ip
        self.running = True

    def run(self):
        print 'Thread Launched'
        first_dir  =  'end-of-1_usrp'+ str(self.ip)
        second_dir =  'end-of-%d_usrp'%proc_time+ str(self.ip)

        self.rx.add_channels(prn_list)
        self.rx.rawfile.seek_rawfile(init_time * 1000 * self.rx.rawfile.S)
        self.rx.scalar_acquisition(prn_list)

        if acq_only:
            self.running = False
            return

        self.rx.scalar_track(mtrack=1000)
        try:
            lock.acquire()
            self.rx.save_measurement_logs(dirname = prepath,subdir= first_dir)
        finally:
            lock.release()

        try:
            self.rx.scalar_track(mtrack=run_time * 1000 - 1000)
            #self.rx.scalar_track(mtrack=39000)
        finally:
            lock.acquire()
            self.rx.save_measurement_logs(dirname = prepath,subdir= second_dir)
            lock.release()
            self.running = False
        self.rx.save_measurement_logs(dirname = prepath,subdir= second_dir)


        self.running = False
        os.makedirs(prepath+'eph%d'%self.ip)
        for prn in self.rx.channels:
            try:
                self.rx.parse_ephemerides(prn_list = [prn],m_start = 40)
                self.rx.channels[prn].ephemerides.save_ephemerides(prepath + 'eph%d/channel%d.mat'%(self.ip,prn))
            except:
                pass
        print 'Scalar Tracking concluded.'
        return

lock= threading.Lock()
s_threads = [scalar_thread(rx,ip_list[i],lock) for i,rx in enumerate(scalar_usrp)]
for st in s_threads:
    st.start()

### Block for scalar
while any([t.running for t in s_threads]):
    print 'Scalar running; total time',run_time
    print 'Current time',[rx._mcount/1000.0 for rx in scalar_usrp]
    time.sleep(30)

print 'Scalar tracking completed.  Launching DP.'

assert (not acq_only)
### Launch DP

dp_usrp = [
    receiver.Receiver(\
        rawfile.RawFile(
            metafile = None,
            abspath  = datpath + prefix[:15] + '_usrp'+str(ip)+'_2500kHz.dat',
            fs = 2.5e6, fi = 0.0e6, ds = 1.0,
            datatype = np.dtype([('i', np.short), ('q', np.short)]),
            notes = 'Data set '+ prefix[:15]
        ), mcount_max = run_time * 50 + 5000
    ) for ip in ip_list
]

for i,dp_rx in enumerate(dp_usrp):
    dp_rx.load_measurement_logs(dirname = prepath, subdir= 'end-of-1_usrp' + str(ip_list[i]))

    del_clist = []
    for prn in dp_rx.channels:
        try:
            dp_rx.channels[prn].ephemerides = scalar_usrp[i].channels[prn].ephemerides
        except:
            del_clist += [prn]
    dp_rx.del_channels(del_clist)

print 'DP Channels'
for i,rx in enumerate(dp_usrp):
    print ip_list[i], rx.channels.keys()

### Time alignment
rxTime_dp_init = []
for rx in dp_usrp:
    rxTime_a, rxTime, posvel_ECEF, posvel_ECI, sats_ECI = naveng.calculate_nav_soln(rx)
    rxTime_dp_init += [rxTime_a]

rxTime_dp_offset = np.round((max(rxTime_dp_init)-np.array(rxTime_dp_init))*1000)
print rxTime_dp_offset
for i,rx in enumerate(dp_usrp):
    rx.scalar_track(mtrack = int(rxTime_dp_offset[i]))

for rx in dp_usrp:
    rx.rawfile.set_rawsnippet_settings(T=0.020,T_big=0.020)
    rx.init_dp()
    print 'Init at',utils.ECEF_to_LLA(rx.ekf.X_ECEF)

### Declare dp threads
keepRunning = True

class rx_thread (threading.Thread):
    def __init__(self,rx, ip, f):
        threading.Thread.__init__(self)
        self.rx = rx
        self.ip = ip

        self.counter = 0
        self.X_list = []
        self.rxTime_list = []
        self.csvfile = f
        self.running = True

    def run(self):
        print 'USRP #',self.ip,'DP Thread Launched'

        printer.header(self.csvfile)
        for mc in range (int((run_time /self.rx.rawfile.T_big))):
            if not keepRunning:
                break
            self.counter += 1
            self.rx.dp_track(1)
            printer.printer(\
                self.counter,\
                weekno,\
                self.rx.rxTime_a,\
                self.rx.ekf.X_ECEF,\
                self.csvfile\
            )

            self.X_list += [self.rx.ekf.X_ECEF.copy()]
            self.rxTime_list += [self.rx.rxTime_a]

            if self.counter % 100 == 0:
                np.save(postpath+'usrp%d_X'%self.ip,np.array(self.X_list))
                np.save(postpath+'usrp%d_t'%self.ip,np.array(self.rxTime_list))
                self.rx.save_measurement_logs(dirname = postpath,subdir= 'end-of-dp_usrp%d'%self.ip)
                print 'DP File saved, continue running.'

        print 'DP Concluded.'
        elapse = time.time() - start
        print elapse,'seconds elapsed for %ds data proc.'%np.ceil(self.rx.rawfile.T_big * mc)
        np.save(postpath+'usrp%d_X'%self.ip,np.array(self.X_list))
        np.save(postpath+'usrp%d_t'%self.ip,np.array(self.rxTime_list))
        #self.rx.save_measurement
        self.rx.save_measurement_logs(dirname = postpath,subdir= 'end-of-dp_usrp%d'%self.ip)
        self.csvfile.close()
        self.running = False

dp_thread = [rx_thread(\
    rx,\
    ip_list[i],\
    open(postpath+'usrp%d.csv'%ip_list[i],'w')\
) for i,rx in enumerate(dp_usrp)]

start = time.time()
for t in dp_thread:
    t.start()

while any([t.running for t in dp_thread]):
    print 'DP running; total time',run_time
    print 'Current time',[t.counter/50.0 for t in dp_thread]
    time.sleep(30)


print 'DP success!'


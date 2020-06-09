
execfile('setting.py')

### Main code starts

from pythonreceiver.libgnss.constants import *
from pythonreceiver.libgnss import rawfile,utils,satpos,ephemeris
from pythonreceiver.scalar import channel, correlator, discriminator, loopfilter,naveng
from pythonreceiver import receiver
import printer
import threading,os

import numpy as np
import scipy.io as sio
import time

import csv


### Launch DP

dp_usrp = [
    receiver.Receiver(\
        rawfile.RawFile(
            metafile = None,
            abspath  = datpath + refname + prefix[:15] + '_usrp'+str(ip)+'_%dkHz.dat'%(fs/1e3),
            fs = fs, fi = 0.0e6, ds = 1.0,
            datatype = np.dtype([('i', np.short), ('q', np.short)]),
            notes = 'Data set '+ refname + prefix[:15]
        ), mcount_max = run_time * 50 + 5000
    ) for ip in ip_list
]

for i,dp_rx in enumerate(dp_usrp):
    dp_rx.load_measurement_logs(dirname = prepath, subdir= 'end-of-1_usrp' + str(ip_list[i]))

    del_clist = []
    for prn in dp_rx.channels:
        try:
            dp_rx.channels[prn].ephemerides = ephemeris.Ephemerides(None)
            dp_rx.channels[prn].ephemerides.load_ephemerides(prepath+'eph%d/channel%d.mat'%(ip_list[i],prn))
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
    rx.rawfile.set_rawsnippet_settings(T=0.020,T_big=0.02)
    rx.init_dp()

    # Load the initialization data from a file to match CUDARecv (overwriting imported logs)
    rx.load_cudarecv_handoff(datpath + cudarecv_handoff + '.csv')
    #rx.perturb_init_ENU(np.matrix([100, 100, 100]).T, rx.ekf.X_ECEF)
    #rx.perturb_init_ECEF(np.matrix([95.2676, 107.907, 3.5931]).T, rx.ekf.X_ECEF)

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
        self.rc_list = []
        self.ri_list = []
        self.fc_list = []
        self.fi_list = []
        self.cp_list = []
        self.rxTime_list = []
        self.csvfile = f
        self.running = True

    def run(self):
        print 'USRP #',self.ip,'DP Thread Launched'

        xFileFile = open(postpath + 'usrp%d_XFile.csv' % self.ip, 'w')
        xwriter = csv.writer(xFileFile, delimiter=',', quotechar='"',
                             quoting=csv.QUOTE_MINIMAL)

        rcFileFile = open(postpath + 'usrp%d_rc.csv' % self.ip, 'w')
        rcwriter = csv.writer(rcFileFile, delimiter=',', quotechar='"',
                             quoting=csv.QUOTE_MINIMAL)
        riFileFile = open(postpath + 'usrp%d_ri.csv' % self.ip, 'w')
        riwriter = csv.writer(riFileFile, delimiter=',', quotechar='"',
                             quoting=csv.QUOTE_MINIMAL)
        fcFileFile = open(postpath + 'usrp%d_fc.csv' % self.ip, 'w')
        fcwriter = csv.writer(fcFileFile, delimiter=',', quotechar='"',
                             quoting=csv.QUOTE_MINIMAL)
        fiFileFile = open(postpath + 'usrp%d_fi.csv' % self.ip, 'w')
        fiwriter = csv.writer(fiFileFile, delimiter=',', quotechar='"',
                             quoting=csv.QUOTE_MINIMAL)
        cpFileFile = open(postpath + 'usrp%d_cp.csv' % self.ip, 'w')
        cpwriter = csv.writer(cpFileFile, delimiter=',', quotechar='"',
                             quoting=csv.QUOTE_MINIMAL)


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
            #temp_rc = []
            #temp_ri = []
            #temp_fc = []
            #temp_fi = []
            #temp_cp = []
            #for prn in self.rx.channels:
            #    temp_rc.append(self.rx.channels[prn].rc[self.rx._mcount])
            #    temp_ri.append(self.rx.channels[prn].ri[self.rx._mcount])
            #    temp_fc.append(self.rx.channels[prn].fc[self.rx._mcount])
            #    temp_fi.append(self.rx.channels[prn].fi[self.rx._mcount])
            #    temp_cp.append(self.rx.channels[prn].cp[self.rx._mcount])
            #self.rc_list += [temp_rc.copy()]
            #self.ri_list += [temp_ri.copy()]
            #self.fc_list += [temp_fc.copy()]
            #self.fi_list += [temp_fi.copy()]
            #self.cp_list += [temp_cp.copy()]

            #self.rxTime_list += [self.rx.rxTime_a]

            xwriter.writerow(np.squeeze(self.rx.ekf.X_ECEF).tolist()[0])
            #rcwriter.writerow(temp_rc)
            #riwriter.writerow(temp_ri)
            #fcwriter.writerow(temp_fc)
            #fiwriter.writerow(temp_fi)
            #cpwriter.writerow(temp_cp)


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
        xFileFile.close()
        rcFileFile.close()
        riFileFile.close()
        fcFileFile.close()
        fiFileFile.close()
        cpFileFile.close()

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


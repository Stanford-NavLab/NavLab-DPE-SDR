#YYYYMMDD_HHMMSS_usrp%d_xxxxkHz.dat

refname    = 'static_opensky_'  # Simulated
datname    = '20180704_190000'
descriptor = 'Talbot_rooftop_4hr_T_big=1'
fs = 2.5e6

cudarecv_handoff = 'handoff_params_usrp6' # Simulated

ip_list       = [6] # Simulated
start_time    = 5 
#  this will run from start_time to (start_time+proc_time)
proc_time     = 30
max_lead_time = 0
weekno        = 2008 # If running simulated data
acq_only      = False
prn_list = [2, 3, 6, 12, 17, 19, 24, 28] # Simulated


datpath  = '/home/ubuntu/repo/'
predir   = './pre-simulator/'
postdir  = './post-simulator/'



### No need to change: house keeping ###

descriptor = descriptor + '_%dUSRP' %len(ip_list)

import os
dir_req  = [datpath,predir,postdir]#,'netcsv','html']
for d in dir_req:
    if not os.path.exists(d):
        os.makedirs(d)

### More housekeeping...

if start_time < 4:
    init_time = start_time
elif start_time < (max_lead_time + 4):
    init_time = 4
else:
    init_time = start_time - max_lead_time

run_time = proc_time  + (start_time - init_time)
prefix   = datname    + '_skip%d_start%d/'%(init_time,start_time)+'6605_7200/'
prepath  = predir     + prefix
postfix  = descriptor + 'proc%ds_6605_7200/'%run_time
postpath = postdir    + postfix
try:
    os.makedirs(postpath)
except:
    if os.listdir(postpath):
        print 'Warning:',postpath,'not an empty directory.'


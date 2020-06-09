from datetime import date
from pythonreceiver import libgnss

def parse_rinex(fname,prn):
    fo  = open(fname, "r")
    #buf = fo.readline()
    while(True):
        buf = fo.readline()
        try:
            p = int(buf[0:2])
        except ValueError:
            pass
        else:
            if (p == prn):
                break
    raw = [buf]
    for i in range(1,8):
        raw += [fo.readline()]
    fo.close()
    eph = libgnss.ephemeris.Ephemerides([])

    val = []
    for i in range(8):
        for j in range(3,79,19):
            if (i == 0 and j < 22):
                continue
            val += [float(raw[i][j:j+19].replace('D','E'))]

    def getTOE(date_o,h,m,s):
        return (date.weekday(date_o)+1)%7*86400+h*3600+m*60+s

    def ORBIT(i,j):
        return val[i*4+j-1]

    eph.prn        = prn

    eph.t_oc       = getTOE(date(int(buf[3:5])+2000,int(buf[6:8]),int(buf[9:11])),\
                     int(buf[12:14]),int(buf[15:17]),float(buf[18:21]))

    eph.a_f0       = ORBIT(0,1);
    eph.a_f1       = ORBIT(0,2);
    eph.a_f2       = ORBIT(0,3);
    eph.IODE       = int(ORBIT(1,0));
    eph.C_rs       = ORBIT(1,1);
    eph.delta_n    = ORBIT(1,2);
    eph.M_0        = ORBIT(1,3);
    eph.C_uc       = ORBIT(2,0);
    eph.e          = ORBIT(2,1);
    eph.C_us       = ORBIT(2,2);
    eph.sqrt_A     = ORBIT(2,3);
    eph.t_oe       = int(ORBIT(3,0));
    eph.C_ic       = ORBIT(3,1);
    eph.OMEGA_0    = ORBIT(3,2);
    eph.C_is       = ORBIT(3,3);
    eph.i_0        = ORBIT(4,0);
    eph.C_rc       = ORBIT(4,1);
    eph.omega      = ORBIT(4,2);
    eph.OMEGADOT   = ORBIT(4,3);
    eph.IDOT       = ORBIT(5,0);
    eph.weeknumber = int(ORBIT(5,2));
    eph.accuracy   = int(ORBIT(6,0));
    eph.health     = int(ORBIT(6,1));
    eph.T_GD       = ORBIT(6,2);
    eph.IODC       = int(ORBIT(6,3));

    eph.timestamp  = {'TOW':506100, 'cp':5190};
    return eph

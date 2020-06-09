# -*- coding: utf-8 -*-


from .constants import *
import numpy as np
import scipy.io as sio
import csv

PARITY_MAT = np.array([[1,1,1,0,1,1,0,0,0,1,1,1,1,1,0,0,1,1,0,1,0,0,1,0],
                       [0,1,1,1,0,1,1,0,0,0,1,1,1,1,1,0,0,1,1,0,1,0,0,1],
                       [1,0,1,1,1,0,1,1,0,0,0,1,1,1,1,1,0,0,1,1,0,1,0,0],
                       [0,1,0,1,1,1,0,1,1,0,0,0,1,1,1,1,1,0,0,1,1,0,1,0],
                       [1,0,1,0,1,1,1,0,1,1,0,0,0,1,1,1,1,1,0,0,1,1,0,1],
                       [0,0,1,0,1,1,0,1,1,1,1,0,1,0,1,0,0,0,1,0,0,1,1,1]])

class Word():
    """
    Instances of the libgnss.Word class are containers for holding
    thirty bits of navigation data, plus the d29 and d30 bits of the previous word.
    """

    def __init__(self, polarity, d29, d30, bits):
        """
        Constructs the libgnss.Word object.

        @type  polarity : int
        @param polarity : This is the polarity of the incoming navigation data bits.
                          It can be either a 1 or a -1. Flip the bits if -1.
                          Note: Polarity is determined by -d30 of the previous word.
                          Note: Polarity of first word is determined by preamble.
        @type  d29 : int
        @param d29 : This is the 29th navigation bit from the previous word.
        @type  d30 : int
        @param d30 : This is the 30th navigation bit from the previous word.
        @type  bits: np.ndarray
        @param bits: This is the np.ndarray containing 30 navigation bits.
        """

        self.polarity = polarity
        self.bits = bits
        self.d29 = bits[28]
        self.d30 = bits[29]
        self.dStarArray = np.array([d29,d30,d29,d30,d30,d29])

        self.bitstring = None
        self.paritypass = None
        self._check_parity()

    def _check_parity(self):
        """
        Helper function that checks the parity of this word.
        """

        assert self.polarity == self.dStarArray[1]
        p = self.dStarArray[1]*PARITY_MAT*self.bits[0:24]
        parities = np.product(np.where(p == 0,1,p),1) * self.dStarArray
        self.paritypass = np.all(parities == self.bits[24:30])
        if self.paritypass:
            bits = np.where(self.polarity*self.bits == -1, 1, 0)
            self.bitstring = ''.join(map('{0:b}'.format, bits))

# Common scaling factors for parsing Subframes

_2_p4 = 2**4
_2_n5 = 2**(-5)
_2_n19 = 2**(-19)
_2_n29 = 2**(-29)
_2_n31 = 2**(-31)
_2_n33 = 2**(-33)
_2_n43 = 2**(-43)
_2_n55 = 2**(-55)

class Subframe():
    """
    Instances of the libgnss.Subframe class are containers for holding
    ten words of navigation data.

    The first word is the TLM which contains the preamble,
    the second word is the HOW which contains the z-count.
    """

    # 'codeperiod' corresponds to the index of the codeperiod at the start of the subframe.
    # 'timeOfWeek' corresponds to the GPS time of week (s) at the start of the subframe,
    #              before the transitition to the first codeperiod.

    def __init__(self, cp, words):
        """
        Constructs a libgnss.Subframe object

        @type  cpcount : int
        @param cpcount : This is the index of the codeperiod that forms the
                         start of the first navbit of the subframe.
        @type  word: list
        @param word: list of libgnss.Word objects
        """

        self.cp = cp
        self.TOW = None

        self.IODE = None
        self.IODC = None
        self.edict = {}
        self.words = words

        self.ID = int(self.words[1].bitstring[19:22], 2)
        self._decode_ephemerides()

    def _decode_ephemerides(self):
        """
        Helper function that decodes the ephemerides contained in this Subframe.

        For more information, see
         - IS-GPS-200H, page 88-120
        """

        if self.ID == 1:
            # IODE and IODC.
            self.IODE = int(self.words[7].bitstring[0:8],2)
            self.IODC = int(self.words[2].bitstring[22:24]+self.words[7].bitstring[0:8],2)
            # TOW
            self.TOW = int(self.words[1].bitstring[0:17], 2)*6 - 6
            self.edict['timestamp'] = {'TOW':self.TOW, 'cp':self.cp}
            print('Subframe'+str(self.ID)+', TOW: %d, cp %d, IODE: %d, IODC: %d'
                  %(self.TOW,self.cp,self.IODE,self.IODC))

            # Ephemeris
            self.edict['weeknumber'] = int(self.words[2].bitstring[0:10],2) + 1024
            self.edict['accuracy'] = int(self.words[2].bitstring[12:16],2)
            self.edict['health'] = int(self.words[2].bitstring[16],2)
            self.edict['T_GD'] = self._twos_complement(self.words[6].bitstring[16:24]) * _2_n31
            self.edict['t_oc'] = int(self.words[7].bitstring[8:24], 2) * _2_p4
            self.edict['a_f2'] = self._twos_complement(self.words[8].bitstring[0:8]) * _2_n55
            self.edict['a_f1'] = self._twos_complement(self.words[8].bitstring[8:24]) * _2_n43
            self.edict['a_f0'] = self._twos_complement(self.words[9].bitstring[0:22]) * _2_n31

        elif self.ID == 2:
            # IODE.
            self.IODE = int(self.words[2].bitstring[0:8], 2)

            # TOW
            self.TOW = int(self.words[1].bitstring[0:17],2)*6 - 6
            self.edict['timestamp'] = {'TOW':self.TOW, 'cp':self.cp}
            print('Subframe'+str(self.ID)+', TOW: %d, cp %d, IODE: %d' %(self.TOW,self.cp,self.IODE))

            # Ephemeris
            self.edict['C_rs'] = self._twos_complement(self.words[2].bitstring[8:24]) * _2_n5
            self.edict['delta_n'] = self._twos_complement(self.words[3].bitstring[0:16]) * _2_n43 * PI
            self.edict['M_0'] = self._twos_complement(self.words[3].bitstring[16:24]
                                                      + self.words[4].bitstring[0:24]) * _2_n31 * PI
            self.edict['C_uc'] = self._twos_complement(self.words[5].bitstring[0:16]) * _2_n29
            self.edict['e'] = int(self.words[5].bitstring[16:24] + self.words[6].bitstring[0:24], 2) * _2_n33
            self.edict['C_us'] = self._twos_complement(self.words[7].bitstring[0:16]) * _2_n29
            self.edict['sqrt_A'] = int(self.words[7].bitstring[16:24] + self.words[8].bitstring[0:24], 2) * _2_n19
            self.edict['t_oe'] = int(self.words[9].bitstring[0:16], 2) * _2_p4

        elif self.ID == 3:
            # IODE
            self.IODE = int(self.words[9].bitstring[0:8], 2)

            # TOW
            self.TOW = int(self.words[1].bitstring[0:17],2)*6 - 6
            self.edict['timestamp'] = {'TOW':self.TOW, 'cp':self.cp}
            print('Subframe'+str(self.ID)+', TOW: %d, cp %d, IODE: %d' %(self.TOW,self.cp,self.IODE))

            # Ephemeris
            self.edict['IDOT'] = self._twos_complement(self.words[9].bitstring[8:22]) * _2_n43 * PI
            self.edict['C_ic'] = self._twos_complement(self.words[2].bitstring[0:16]) * _2_n29
            self.edict['OMEGA_0'] = self._twos_complement(self.words[2].bitstring[16:24] + self.words[3].bitstring[0:24]) \
                                                          * _2_n31 * PI
            self.edict['C_is'] = self._twos_complement(self.words[4].bitstring[0:16]) * _2_n29
            self.edict['i_0'] = self._twos_complement(self.words[4].bitstring[16:24] + self.words[5].bitstring[0:24]) \
                                                      * _2_n31 * PI
            self.edict['C_rc'] = self._twos_complement(self.words[6].bitstring[0:16]) * _2_n5
            self.edict['omega'] = self._twos_complement(self.words[6].bitstring[16:24] + self.words[7].bitstring[0:24]) \
                                                        * _2_n31 * PI
            self.edict['OMEGADOT'] = self._twos_complement(self.words[8].bitstring[0:24]) * _2_n43 * PI

        elif self.ID == 4:
            # TOW
            self.TOW = int(self.words[1].bitstring[0:17], 2)*6 - 6
            self.edict['timestamp'] = {'TOW':self.TOW, 'cp':self.cp}
            print('Subframe'+str(self.ID)+', TOW: %d, cp %d' %(self.TOW,self.cp))

        elif self.ID == 5:
            # TOW
            self.TOW = int(self.words[1].bitstring[0:17], 2)*6 - 6
            self.edict['timestamp'] = {'TOW':self.TOW, 'cp':self.cp}
            print('Subframe'+str(self.ID)+', TOW: %d, cp %d' %(self.TOW,self.cp))

        else:
            print('Ephemerides decoding unsuccessful: Unknown subframe ID: '+str(self.ID))

    def _twos_complement(self, s):
        """
        Returns the binary-to-decimal two's-complement of the binary string,
        see http://en.wikipedia.org/wiki/Two's_complement.

        @type  s : str
        @param s : This is the binary string to convert to decimal.
        @rtype  : int
        @return : This is the integer value conversion of the binary string.
        """

        decimal = int(s, 2) if s[0] == '0' else int(s, 2) - 2**len(s)
        return decimal


# Ephemerides attributes

_CLOCK_NAMES = ['weeknumber', 'accuracy', 'health', 'T_GD', 't_oc', 'a_f2',
                'a_f1', 'a_f0']

_EPHEMERIS_NAMES = ['C_rs', 'delta_n', 'M_0', 'C_uc', 'e', 'C_us',
                    'sqrt_A', 't_oe', 'C_ic', 'OMEGA_0', 'C_is', 'i_0',
                    'C_rc', 'omega', 'OMEGADOT', 'IDOT']

_NONSTANDARD_NAMES = ['timestamp','total','complete','IODE','IODC']

_NAMES = _CLOCK_NAMES + _EPHEMERIS_NAMES + _NONSTANDARD_NAMES
_TOTAL = len(_NAMES)

class Ephemerides():
    """
    Instances of the libgnss.Ephemerides class store ephemerides
    unique to a particular issue (IODE and IODC).
    """

    def __init__(self, subframes):
        """
        Constructs a libgnss.Ephemerides object.

        For more information, see
         - IS-GPS-200H, page 88-120

        @type  subframes : list
        @param subframes : list of libgnss.Subframe objects
        """
        self.subframes = subframes

        # Set up the issues
        self.IODE = None
        self.IODC = None

        # Set up the ephemerides attributes
        for name in _NAMES:
            setattr(self, name, None)

        # Create a flag to let us know how many ephemerides are set.
        self.total = 0
        self.complete = False

        if subframes is None:
            return

        for subframe in self.subframes:
            if subframe.ID in [1,2,3]:
                if self.IODE is None:
                    self.IODE = subframe.IODE
                if (self.IODC is None) and (subframe.ID == 1):
                    self.IODC = subframe.IODC
                if self.IODE == subframe.IODE:
                    self.add(subframe.edict)
                else:
                    print('Ephemerides decoding unsuccesful: IODE has changed.')

        self.complete = self.is_complete()

    def add(self, edict):
        """
        Add new ephemeris parameters specified in edict.

        @type  edict: dict
        @param edict: dictionary of ephemeris parameters.
        """

        # Fill in any ephemerides that are passed in kwargs.
        for name in edict.keys():
            assert name in _NAMES
            if (self.__dict__[name] is None) and (edict[name] is not None):
                self.total = self.total + 1
                self.__dict__[name] = edict[name]
            else:
                if (self.__dict__[name] is not None) and (edict[name] is not None):
                    print('Warning: Ephemerides already contains ephemerides parameter: '+name)
                elif (edict[name] is None):
                    print('Warning: Subframe does not contain ephemerides parameter: '+name)
                else:
                    print('Warning: Unknown ephemerides addition error.')

    def is_complete(self):
        """
        Return True if the Ephemerides object is full (ready for use).
        @rtype: bool
        @return: Return True if the Ephemerides object is full (ready for use).
                 Return False if the Epehemerides object is not full.
        """
        return all([self.IODE!= None, self.IODC != None, self.total == _TOTAL])

    def print_ephemerides(self):
        """
        Prints the ephemerides
        """
        print('IODE: '+str(self.IODE)+', IODC: '+str(self.IODC))
        for name in _NAMES:
            print(name+': '+str(self.__dict__[name]))

    def save_ephemerides(self, filename, csvfname):
        save_dict= {}
        for name in _NAMES:
            save_dict[name] = getattr(self,name)
        sio.savemat(filename,save_dict)
        with open(csvfname, 'wb') as csv_file:
            writer = csv.writer(csv_file)
            for key, value in save_dict.items():
                writer.writerow([key, value])


    def save_scalar_handoff_ephem(self, prn_list, dirname=None, subdir=''):

        save_dict= {}
        for name in _NAMES:
            temp = getattr(self,name)
            # Iterate through the list to order params according to order of prn_list
            for prn in prn_list:
                save_dict[name].append(temp[prn])

        with open(os.path.join(dirname, subdir, 'handoff_params.csv'), 'a') as csv_file:
            writer = csv.writer(csv_file)
            for key, value in save_dict.items():
                temp = [key]
                temp.extend(value)
                writer.writerow(temp)


            print('Saved ephems params')
        return

    def load_ephemerides(self, filename):
        load_dict = sio.loadmat(filename)
        for name in _NAMES:
            setattr(self,name,load_dict[name][0,0])

        # Loading in ephemerides in the file provided loads in a variable named timestamp.
        # This variable is picked apart in the following lines to extract 'cp' and 'TOW'
        self.timestamp = {\
                'cp':self.timestamp[0][0,0],\
                'TOW':self.timestamp[1][0,0]\
        }

        self.complete = True if self.complete else False

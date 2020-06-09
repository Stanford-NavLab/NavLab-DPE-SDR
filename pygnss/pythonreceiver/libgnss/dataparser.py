# -*- coding: utf-8 -*-


from . import ephemeris as eph
import numpy as np

preamble = np.array([-1,1,1,1,-1,1,-1,-1])
preamble_cp = np.kron(preamble, np.ones(20))
        
def parse_ephemerides(channel, m_start, m_end):
    
    cp_start = int(channel.cp[m_start])
    cp_end = int(channel.cp[m_end])
    
    if len(np.where(np.diff(channel.cp[m_start:m_end])!=1)[0]) != 0:
        print('WARNING: possible code period slip occured.')

    cp = np.arange(cp_start,cp_end)
    iP = channel.cp_sign[cp]

    preamble_correlations = np.correlate(iP,preamble_cp,'valid')
    preamble_locations = np.where(np.abs(preamble_correlations)>153)[0][:]
    preamble_locations_set = set(preamble_locations)
    subframe_found = False

    for test in preamble_locations:
        test_set = set([test,test+6000,test+6000*2,test+6000*3,test+6000*4])
        print(test_set,test_set.issubset(preamble_locations_set))
        if test_set.issubset(preamble_locations_set):
            subframe_found = True
            subframe_locations = sorted(list(test_set))
            break

    if subframe_found == False:
        print('Ephemerides decoding unsuccessful: Unable to locate 5 preambles within given range.')
        return

    subframe_polarities = np.sign(preamble_correlations[subframe_locations])
    print('polarities: '+str(subframe_polarities))

    if not (np.any(subframe_polarities == [-1,-1,-1,-1,-1]) or np.any(subframe_polarities == [1,1,1,1,1])):
        print('Ephemerides decoding warning: Bit flip occured in between subframes')

    navbits_all = np.reshape(iP[subframe_locations[0]:(subframe_locations[0]+6000*5)],(1500,20))
    navbits_all = np.sign(np.sum(navbits_all,1))

    subframes = []
    first_d29d30 = np.reshape(iP[(subframe_locations[0]-40):(subframe_locations[0])],(2,20))
    first_d29d30 = np.sign(np.sum(first_d29d30,1))
    d29 = first_d29d30[0]
    d30 = first_d29d30[1]

    reload(eph)
    for subframeNr in range(5):
        words = []
        wordNr = 0    
        polarity = subframe_polarities[subframeNr]
        bits = navbits_all[(subframeNr*300+(wordNr*30)):(subframeNr*300+(wordNr*30+30))]     
        words.append(eph.Word(polarity,d29,d30,bits))
        for wordNr in range(1,10):
            bits = navbits_all[(subframeNr*300+(wordNr*30)):(subframeNr*300+(wordNr*30+30))]
            words.append(eph.Word(words[-1].d30,words[-1].d29,words[-1].d30,bits))
        d29 = words[-1].d29
        d30 = words[-1].d30
        subframes.append(eph.Subframe(subframe_locations[subframeNr]+cp[0],words))

    parsed_ephemerides  = eph.Ephemerides(subframes)
    channel.ephemerides = parsed_ephemerides
    
    return

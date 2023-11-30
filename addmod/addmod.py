import os
import builtins      as bt
import numpy         as np
import basics.basics as bs
import bcfits.bcfits as bc
import bcfilt.bcfilt as fl
import bchead.bchead as hd
import common.common as cn

def pyadd_bursts(file1):
    # python version of add_bursts code
    global iread, b1, a1, b2, a2, bp1, bp2
    global x, t, y1,t1,m1,p1,s1, y2,t2,m2,p2,s2

    #Init variables
    bt.iread = True	# Read filter file only on first call

    # Ask for burst parameters and store for later usage
    if not bt.jupy:
        b1,a1,f1,b2,a2,f2,fo = burstparams(file1)
    else:
        b1,a1,f1,b2,a2,f2,fo = file1
    bt.bp1 = [b1,a1,f1]
    bt.bp2 = [b2,a2,f2]

    # Read BC/CB model in first fits file
   #x,y1,t1,e1,s1,m1,p1,*nouse = bc.read_bcfits(f1)
    x,y1,t1,e1,s1,m1,p1,*nouse = bc.bcfits(f1)
    hd.multhead(0,f1,t1,x,e1)

    # Read BC/CB model in second fits file
    if f2 != f1:
       #w2,y2,t2,e2,s2,m2,p2,*nouse = bc.read_bcfits(f2)
        w2,y2,t2,e2,s2,m2,p2,*nouse = bc.bcfits(f2)
        if len(w2) != len(x):
            print(' The two models must use the same spectral library.')
            quit()
        hd.multhead(0,f2,t2,w2,e2)
    else:
        y2,t2,m2,p2,s2 = y1,t1,m1,p1,s1
    bt.s1 = s1
    bt.s2 = s2

    # Build combined time scale for the two bursts
    t=[]
    for i in range(len(t1)):
        if t1[i] <= 14.E9:
            t.append(t1[i])
    for i in range(len(t2)):
        if b2+t2[i] <= 14.E9:
            t.append(t2[i]+b2)
    t = np.unique(t)            # unique sorts array t and suppresses duplicate entries (no need for t = np.sort(np.unique(t)))
    t = checkunique(t,1.E-6)    # suppress time steps differing by less than eps = 1.E-6 in the log

    # Add 2 bursts
    addbursts('ADD_BURSTS',fo)

def burstparams(file1):
    # Ask for input parameters for 2 bursts
    print()
    print (' Galaxy Spectral Evolution Library (GALAXEV)')
    print (' python version (C) 2019-2023 - G. Bruzual and S. Charlot - All Rights Reserved')
    print()
    print('  This program computes the combined sed for 2 bursts of star formation.')
    print('  The sed corresponding to each burst must be computed first. Each burst')
    print('  can follow an arbitrary star formation law. A burst is characterized by')
    print('  the beginning time (in Gyr) and the burst strength. The FIRST burst is')
    print('  assumed to start at t = 0 with strehgth = 1. The strength of the SECOND')
    print('  burst is then relative to that of the FIRST burst.')
    print()
    f1 = input(' First BC_GALAXEV sed in file [' + bs.lnam(file1) + ']  = ')
    if len(f1) <= 0:
        f1 = bs.lnam(file1)
    else:
        f1 = bs.zrep(file1,f1)
    f1 = bs.fcheck(f1)
    f2 = input(' Second BC_GALAXEV sed in file [' + bs.lnam(f1) + '] = ')
    if len(f2) <= 0:
        f2 = bs.lnam(f1)
    else:
        f2 = bs.zrep(f1,f2)
    f2 = bs.fcheck(f2)
    print()
   #b1 = input(' Burst 1: Beginning time (Gyr), Burst amplitude = ')
   #b1 = b1.replace(',',' ') ; b1= b1.split() ; a1 = float(b1[1]) ; b1 = float(b1[0])*1.E9
    b1 = 0. ; a1 = 1.
    print(' Burst 1: Beginning time (Gyr), Burst amplitude = 0,1')
    b2 = input(' Burst 2: Beginning time (Gyr), Burst amplitude = ')
    b2 = b2.replace(',',' ') ; b2= b2.split() ; a2 = float(b2[1]) ; b2 = float(b2[0])*1.E9
    print()
    fo = input(' Output file name = ')
    if fo=='q':
        quit()
    return b1,a1,f1,b2,a2,f2,fo

def checkunique(t,eps):
    # Keeps only time steps differing in more than eps
    tx=[]
    tx.append(t[0])
    tx.append(t[1])
    for i in range(2,len(t)):
        tl = np.log10(np.float64(t[i])) - np.log10(np.float64(t[i-1]))
        if tl > eps:
            tx.append(t[i])
    return tx

def addbursts(ly,fo):
    # Add 2 bursts and compute model properties
    global iread, kf, kp, fd

    # Init variables
    p=-1		# Needed by 'percent' function on first call
    ly = ' ' + ly + ' --> ' + fo   # header for percent script

    for i in range(len(t)):
        if bt.iread:
            # Init output tables
            bt.kf,kp = cn.opent(1)
            # Read filter file and build filter arrays at z = 0
            bt.ffd = fl.qfilters(x,kf,0.,0,1)
            # Store model SED wavelength array
            bt.td1.add_column(np.array(x))

        # Add spectra, interpolating when needed
        b = []
        if t[i] < b2:
            y = a1*bc.ft(t[i],y1,t1)
            for j in range(len(kp)):
                b.append(a1*bc.pt(t[i],p1,m1,t1,kp[j]))
            sf = a1*bc.pt(t[i],p1,m1,t1,23)
        else:
            y = a1*bc.fi(t[i],y1,t1) + a2*bc.fi(t[i]-b2,y2,t2)
            for j in range(len(kp)):
                b.append(a1*bc.pi(t[i],p1,m1,t1,kp[j]) + a2*bc.pi(t[i]-b2,p2,m2,t2,kp[j]))
            sf = a1*bc.pi(t[i],p1,m1,t1,23) + a2*bc.pi(t[i]-b2,p2,m2,t2,23)

        # Compute model rest-frame properties at this age
        cn.rf_props(t[i],x,y,b,sf,1)

        # Store model SED = y
        bt.td1.add_column(np.array(y))

        # Report percent done
        p = cn.percent(i,len(t),ly,p)

    # Write results to fits file
    bc.wrfits(fo,1)


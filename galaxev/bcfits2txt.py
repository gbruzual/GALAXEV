#!/usr/bin/env python3

# G. Bruzual (June 2022)

import os
import sys
import math
import numpy as np
from   scipy.io import FortranFile

def lnam(s):
    # Returns file name without directory path
    if '/' in s:
        s = s.split('/')
        return s[len(s)-1]
    else:
        return s

def myheader(f,e,l):
    # Writes header in BC tables
    a = '#          '
    h = ['3', '3', '3', '3', '5', '5', '5', '3', '3', '3']
    i = 0
    for k in range(1000):
        if 'END' in e['HEADER' + str(i+1)]:
           break
        else:
            s = e['HEADER' + str(i+1)] + e['HEADER' + str(i+2)]
            s = s.replace('|',' ')
            if k==1:
                u = len(s)
            if k==3 and len(s) > u:
                s = s.replace('  I','I')
            f.write(a + s + '\n')
            i=i+2 
    f.write('#' + '\n')
    if l == 0:
        # Header for file .0VEGAmag to write Vega magnitudes
        h1 = '#    (1)       (2)       (3)       (4)       (5)       (6)       (7)       (8)       (9)       (10)      (11)      (12)      (13)      (14)      (15)      (16)      (17)'
        h2 = '#'
        h3 = "#log-age-yr    Mbol       U         B2        B3        V         Rc        Ic        R         I         J         K         L        PalJ      PalH      PalK       K'"

    elif l == 1:
        # Header for file .1VEGAmag to write Vega magnitudes
        h1 = '#    (1)       (2)       (3)       (4)       (5)       (6)       (7)       (8)       (9)       (10)      (11)      (12)      (13)      (14)      (15)      (16)      (17)'
        h2 = '#'
        h3 = '#log-age-yr    Mbol      Vmag     2MASSJ    2MASSH   2MASSKs   IRAC3.5   IRAC4.5   IRAC5.7   IRAC7.9    IRAS12    IRAS25    IRAS60   IRAS100    MIPS24    MIPS70   MIPS160'

    elif l == 2:
        # Header for file .2ABmag to write AB mag (SDSS, CFHT MegaCam, and GALEX)
        h1 = '#    (1)       (2)       (3)       (4)       (5)       (6)          (7)       (8)       (9)       (10)      (11)      (12)      (13)      (14)      (15)      (16)         (17)      (18)         (19)      (20)      (21)         (22)         (23)'
        h2 = '#             <---------------- SDSS AB mag --------------->       <---------------------- CFHT MegaCam 1st and 3rd generation filters AB mag -------------------->       2Mass    CFHT Ks      <---------------- GALEX AB mag and flux --------------->'
        h3 = '#log-age-yr     u         g         r         i         z            u1        u3        g1        g3        r1        r3        i2        i3        z1        z3           H         Ks          FUV       NUV      F(FUV)       F(NUV)       F(1500A)'

    elif l == 3:
        # Header for file .3ABmag to write AB mag (HST WFC and ACS bands)
        h1 = '#    (1)       (2)       (3)       (4)       (5)       (6)       (7)        (8)       (9)      (10)      (11)      (12)      (13)      (14)      (15)      (16)      (17)      (18)      (19)      (20)      (21)      (22)      (23)      (24)      (25)      (26)      (27)      (28)      (29)      (30)      (31)      (32)      (33)'
        h2 = '#             <------------------------------------------- HST WFC3 UVIS1 LEGUS ABmag ------------------------------------------->    <----------------------------------- HST WFC3 ABmag --------------------------------->    <------------------------------------------- HST ACS WFC ABmag ------------------------------------------>'
        h3 = '#log-age-yr    F225w     F275w     F336w     F438w     F547m     F555w     F606w     F625w     F656n     F657n     F658n     F814w     F110W     F125W     F160W     F225W     F336W    FR388N     F438W     F555W     F814W     F220w     F250w     F330w     F410w     F435w     F475w     F555w     F606w     F625w     F775w     F814w'

    elif l == 4:
        # Header for file .4lsindx to write Lick indices measured from sed (1-21)
        h1 = '#    (1)       (2)      (3)      (4)      (5)      (6)      (7)      (8)      (9)      (10)     (11)     (12)     (13)     (14)     (15)     (16)     (17)     (18)     (19)     (20)     (21)     (22)'
        h2 = '#'
        h3 = '# Index_No.      1:       2:       3:       4:       5:       6:       7:       8:       9:      10:      11:      12:      13:      14:      15:      16:      17:      18:      19:      20:      21:'
        h4 = '# log-age      CN_1     CN_2   Ca4227    G4300   Fe4383   Ca4455   Fe4531   Fe4668    Hbeta   Fe5015     Mg_1     Mg_2     Mg-b   Fe5270   Fe5335   Fe5406   Fe5709   Fe5782     Na-D    TiO_1    TiO_2'
        h5 = '#   (yr)      (mag)    (mag)     (A)      (A)      (A)      (A)      (A)      (A)      (A)      (A)     (mag)    (mag)     (A)      (A)      (A)      (A)      (A)      (A)      (A)     (mag)    (mag)'

    elif l == 5:
        # Header for file .5lsindx to write other indices measured from sed
        h1 = '#    (1)        (2)      (3)      (4)      (5)      (6)      (7)       (8)        (9)        (10)      (11)      (12)      (13)      (14)      (15)'
        h2 = '#'
        h3 = '# Index_No.    WO-1:    WO-2:    WO-3:    WO-4:    GC-1:    4000A    DTT-Ca1    DTT-Ca2    DTT-Ca3   DTT-MgI     DM-04     DM-04     DM-04       BH:'
        h4 = '# log-age    HdeltaA  HgammaA  HdeltaF  HgammaF  D(4000)    B4_VN   CaII8498   CaII8542   CaII8662   MgI8807   H8_3889   H9_3835  H10_3798     BH-HK'
        h5 = '#   (yr)        (A)      (A)      (A)      (A)       .        .         (A)        (A)        (A)       (A)       (A)       (A)       (A)        (A)'

    elif l == 6:
        # Header for file .6lsindx to write Fanelli et al. UV spectral indices measured from sed
        h1 = '#'
        h2 = '#    (1)        (2)        (3)        (4)        (5)        (6)        (7)        (8)        (9)        (10)       (11)       (12)       (13)       (14)       (15)       (16)       (17)       (18)       (19)       (20)       (21)       (22)'
        h3 = '#'
        h4 = '# log_t(yr)    BL1302      SiIV      BL1425     Fe1453    CIV1548a   CIV1548c   CIV1548e    BL1617     BL1664     BL1719     BL1853    FeII2402    BL2538    FeII2609     MgII       MgI       Mgwide      FeI       BL3096     CIVabs     HeIIems'
        h5 = '#   (yr)        (A)        (A)        (A)        (A)        (A)        (A)        (A)        (A)        (A)        (A)        (A)        (A)        (A)        (A)        (A)        (A)        (A)        (A)        (A)        (A)        (A)'

    elif l == 7:
        # Header for file .7physp to write model properties
        h1 = '#    (1)       (2)     (3)      (4)      (5)        (6)      (7)      (8)      (9)      (10)        (11)        (12)        (13)        (14)        (15)        (16)        (17)        (18)        (19)'
        h2 = '#'
        h3 = '#log-age-yr   B4000    B4_VN  B4_SDSS    B912       NLy      NHeI    NHeII  NHeII/NLy   SNe/yr     PISNe/yr    RegSNe/yr   IaSNe/yr  FailedSNe/yr  PNBR/yr      N(BH)       N(NS)       N(WD)       Lx(Lo)'
    
    elif l == 8:
        # Header for file .8physp to write model properties
        h1 = '#    (1)       (2)       (3)       (4)        (5)       (6)          (7)          (8)          (9)         (10)         (11)           (12)        (13)        (14)        (15)        (16)        (17)'
        h2 = '#                                                                                                                      M*_tot='
        h3 = '#log-age-yr    Mbol      Bmag      Vmag      Kmag      M*_liv     M_remnants    M_ret_gas    M_galaxy      SFR/yr    M*_liv+M_rem   M*_tot/Lb   M*_tot/Lv   M*_tot/Lk   M*_liv/Lb   M*_liv/Lv   M*_liv/Lk'

    elif l == 9:
        # Header for file .9physp to write model properties
        h1 = '#    (1)       (2)         (3)          (4)          (5)          (6)          (7)          (8)'
        h2 = '#                                                                          Total Mass   Total Dust'
        h3 = "#log-age-yr    Mbol    b(t)*'s/yr   B(t)/yr/Lo  Turnoff_mass   BPMS/BMS     Loss Rate   Prod. Rate"

    else:
        np = str(len(t))
        nl = str(len(w))
        h1 = '#          The first record contains the time scale (' + np + ' steps). The second record contains the wavelength scale (' + nl + ' points)'
        h2 = '#          The next ' + np + ' records contain the SED for each time step (' + nl + ' points per SED)'
        h3 = '#          Then 20 records follow with various physical properties vs. age (' + np + ' points per record)'

    f.write(h1 + '\n')
    f.write(h2 + '\n')
    f.write(h3 + '\n')
    if h[l] == '5':
        f.write(h4 + '\n')
        f.write(h5 + '\n')
    if l==-1:
        f.write('#' + '\n')

def dark(a,of):
    # Returns number of time steps in which the galaxy is still dark
    for i in range(len(t)-len(a)):
        if i==0:
            s = '#%9.6f' % 0.00 + '    Dark galaxy'
        else:
            s = '#%9.6f' % math.log10(t[i]) + '    Dark galaxy'
        of.write(s + '\n')
    return

def xtab():
    # Writes header in BC tables
#   a = '#          '
    x = ['.0VEGAmag', '.1VEGAmag', '.2ABmag', '.3ABmag', '.4lsindx', '.5lsindx', '.6lsindx', '.7physp', '.8physp', '.9physp']
    f = lnam(ifile.replace('.fits',''))
    o = '.'
    print('            Creating files: ' + f + '.0VEGAmag + (.1VEGAmag, .2ABmag, .3ABmag, .4lsindx, .5lsindx, .6lsindx, .7physp, .8physp, .9physp)')
    print()
    o = o + '/' + f
    a = m['logage']				# time scale without pre MS evolution (agrees with t for BC03 and C&B models)
    for l in range(10):
        of = o + x[l]
        of = open(of,'w')
        myheader(of,e,l)

        if l == 0:
            # Vega mags
            h = [ 2, 14, 23, 22, 32, 43, 49, 44, 53, 59, 69, 71, 61, 63, 70, 66]
            u = [0]*16
            dark(a,of)
            for r in range(len(a)):
                for j in range(16):
                    u[j] = m[r][h[j]-1]
                s = '%10.6f' % a[r] + ' ' + ' '.join('%9.4f' % v for v in u)
                of.write(s + '\n')
            of.close()

        if l == 1:
            # Vega mags
            h = [ 2, 32, 58, 65, 68, 72, 73, 74, 75, 76, 78, 79, 81, 77, 80, 82]
            u = [0]*16
            dark(a,of)
            for r in range(len(a)):
                for j in range(16):
                    u[j] = m[r][h[j]-1]
                s = '%10.6f' % a[r] + ' ' + ' '.join('%9.4f' % v for v in u)
                of.write(s + '\n')
            of.close()

        if l == 2:
            # AB mags
            h = [13, 24, 35, 45, 55, 16, 15, 27, 26, 36, 39, 46, 47, 54, 56, 64, 67, 3, 5, 83, 84, 85]
            u = [0]*22
            dark(a,of)
            for r in range(len(a)):
                for j in range(22):
                    u[j] = m[r][h[j]-1]
                s = ('%10.6f' % a[r]
                     +    ' ' + ' '.join('%9.4f'  % z for z in u[0:5])
                     + '    ' + ' '.join('%9.4f'  % z for z in u[5:15])
                     + '    ' + ' '.join('%9.4f'  % z for z in u[15:17])
                     + '    ' + ' '.join('%9.4f'  % z for z in u[17:19])
                     + ' '    + ' '.join('%12.4E' % z for z in u[19:22]) )
                of.write(s + '\n')
            of.close()

        if l == 3:
            # AB mags
            h = [ 7, 8, 11, 21, 31, 28, 34, 37, 40, 41, 42, 50, 57, 60, 62, 6, 12, 18, 20, 29, 52, 4, 9, 10, 17, 19, 25, 30, 33, 38, 48, 51]
            u = [0]*32
            dark(a,of)
            for r in range(len(a)):
                for j in range(32):
                    u[j] = m[r][h[j]-1]
                s = '%10.6f' % a[r] + ' ' + ' '.join('%9.4f' % v for v in u)
                of.write(s + '\n')
            of.close()

        if l == 4:
            # Lick indices
            u = [0]*21
            dark(a,of)
            for r in range(len(a)):
                for j in range(1,22):
                    u[j-1] = d[r][j]
                s = '%10.6f' % a[r] + ' ' + ' '.join('%8.4f' % v for v in u)
                of.write(s + '\n')
            of.close()

        if l == 5:
            # Other indices
            u = [0]*14
            dark(a,of)
            for r in range(len(a)):
                for j in range(22,36):
                    u[j-22] = d[r][j]
                s = ('%10.6f' % a[r]
                     + '  ' + ' '.join('%8.4f'  % z for z in u[0:6])
                     +  ' ' + ' '.join('%10.4f' % z for z in u[6:9])
                     +  ' ' + ' '.join('%9.4f'  % z for z in u[9:14]) )
                of.write(s + '\n')
            of.close()

        if l == 6:
            # UV indices
            u = [0]*21
            dark(a,of)
            for r in range(len(a)):
                for j in range(36,57):
                    u[j-36] = d[r][j]
                s = '%10.6f' % a[r] + ' ' + ' '.join('%10.4f' % v for v in u)
                of.write(s + '\n')
            of.close()

        if l == 7:
            # Physical properties
            h = [ 9, 7, 8, 6, 2, 3, 4, 5, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
            u = [0]*18
            dark(a,of)
            for r in range(len(a)):
                for j in range(18):
                    u[j] = p[r][h[j]-1]
                s = (  '%10.6f' % a[r]
                     +    ' ' + ' '.join('%8.4f'  % z for z in u[0:3])
                     + '%11.3E' % u[3]
                     + ' '    + ' '.join('%8.4f'  % z for z in u[4:8])
                     + ' '    + ' '.join('%11.4E' % z for z in u[8:19]) )
                of.write(s + '\n')
            of.close()

        if l == 8:
            # Physical properties
            h = [ 2, 22, 32, 70, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31]
            u = [0]*16
            dark(a,of)
            for r in range(len(a)):
                for j in range(4):
                    u[j] = m[r][h[j]-1]
                for j in range(4,16):
                    u[j] = p[r][h[j]-1]
                s = (  '%10.6f' % a[r]
                     +    ' ' + ' '.join('%9.4f'  % z for z in u[0:4])
                     + ' '    + ' '.join('%12.4E' % z for z in u[4:10])
                     + '  '   + ' '.join('%11.4E' % z for z in u[10:17]) )
                of.write(s + '\n')
            of.close()

        if l == 9:
            # Physical properties
            h = [ 2, 32, 33, 34, 35, 36, 37]
            u = [0]*7
            dark(a,of)
            for r in range(len(a)):
                for j in range(1):
                    u[j] = m[r][h[j]-1]
                for j in range(1,7):
                    u[j] = p[r][h[j]-1]
                s = '%10.6f' % a[r] + '%10.4f' % u[0] + ' ' + ' '.join('%12.4E' % z for z in u[1:8])
                of.write(s + '\n')
            of.close()

def getr(n,k):
    # Builds age array for a range selected by record numbers
    n  = n.replace(',',' ')
    n  = n.replace(':',' ')
    n  = n.split()
    n1 = abs(int(n[0]))
    n2 = abs(int(n[1]))+1
    n = [0]*(n2-n1)
    i = -1
    for j in range(n1,n2):
        i=i+1
        n[i] = t[j-1]*1.E-9
    return n,k

def geti(a,t):
    # Find number of records i[k] corresponding to age a[k] (in Gyr) selected by user
    n = len(a)				# Number of time steps requested
    i = [""]*(n+1)
    k = 0
    for l in range(n):
        k = k+1
        c = a[l]
        if c > 1.E3:
            c = c*1.E-9			# if age entered in yr transform to Gyr
        j = np.abs(t-c*1.E9).argmin()
        i[k] = int(j)+1
    k1 = int(i[1])-1
    k2 = int(i[len(i)-1])-1
    if k1 == 0 and k2 == 220:
        x = 'sed'
    elif k2 == k1:
        x = f'sed.{t[k1]:9.3E}yr'.replace('+','')
    else:
        x = f'sed.{t[k1]:9.3E}_{t[k2]:9.3E}yr'.replace('+','')
    return i,x

def examples():
    # Shows some examples of how to use this procedure
    print ('')
    print ('                  Examples: Screen input                       Command line                                              Output')
    print ("                    Choice: 0.001 0.05 0.1 1 10 11 12 13 14    bcfits2txt file.fits '0.001 0.05 0.1 1 10 11 12 13 14'    SEDs t = 1Myr, 50 Myr, 100Myr, 1, 10, 11, 12, 13, 14 Gyr")
    print ("                            -1 10 20 50 100 150 175 200 220    bcfits2txt file.fits '-1 10 20 50 100 150 175 200 220'    SEDs N = 1, 10, 20, 50, 100, 150, 175, 200, 220")
    print ('                            200:210                            bcfits2txt file.fits 200:210                              SEDs N = 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210')
    print ('                            a                                  bcfits2txt file.fits a                                    SEDs N = 1 to 221')
    print ('                            t                                  bcfits2txt file.fits t                                    Tables with model properties')
    print ('                            i                                  bcfits2txt file.fits i                                    create file.ised2 file')
    print ('                            c                                  bcfits2txt file.fits c                                    create file.ised_ASCII file')
    print ()

def deca(a):
    # Decode records to write in file
    a = a.replace(',',' ')
    if a == 'q' or len(a) <= 0:
        sys.exit()
    elif a == 't':
        k = 0
        return a,k
    elif a == 'i':
        # write ised file
        k = -2
        return a,k
    elif a == 'c':
        # write ised_ASCII file
        k = -3
        return a,k
    elif a == 'h':
        examples()
        k = -1
        return a,k
    elif a == 'a':
        a = '1:' + str(len(t))
    if (a.find(':') != -1):
        # if a=n1:n2, extract sed by record number from n1 to n2
        a = a.replace('-','')
        a,k = getr(a,2)
        return a,k
    if (a[0] == '-'):
        # if first number is < 0, extract sed by record number
        a = a[1:]
        a = a.replace('-','')
        a = a.split()
        a = np.array(a, dtype=np.int32)
        a = np.unique(a)		# unique sorts array a and suppresses duplicate entries (no need for a = np.sort(np.unique(a)))
        b = [0]*len(a)
        for i in range(len(a)):
            b[i] = t[a[i]-1]*1.E-9
        k = 2
        return b,k
    else:
        # Extract by age in Gyr
        a = a.split()
        a = np.array(a, dtype=np.float32)
        for i in range(len(a)):	# if age entered in yr transform to Gyr
            if a[i] > 1.E3:
                a[i] = a[i]*1.E-9
        a = np.unique(a)		# unique sorts array a and suppresses duplicate entries (no need for a = np.sort(np.unique(a)))
        k = 1
        return a,k

def getc(a):
    # Ask for time steps to plot
    a = input(a)
    a,k = deca(a)
    return a,k

def getb(a):
    # Check command line argument
    b = ''.join(str(e) for e in a)
    b,k = deca(b)
    return b,k

def nsed(n):
    # Build file name and open output file
    i,x = geti(n,t)
    o   = ifile.replace('fits',x)
    of  = open(o,'w')
    print('\n               Writing records:'  + ''.join('%4i' % i[v] for v in range(1,len(i))))
    print('                       to file: ' + o + '    (will take some time...)')
    # Write output file header
    b = [0]*len(i)
    s = ['']*len(i)
    for j in range(1,len(i)):
       #b[j] = a[int(i[j])-1]
        b[j] = t[int(i[j])-1]
        s[j] = ' at age yr'
    of.write('# Input file name  = ' + ifile + '\n')
    of.write('# Output file name = ' + o     + '\n')
    of.write('# Record'     + ' '.join('%10i'       % i[v] for v in range(1,len(i)))   + '\n')
    of.write('# Column'     + ' '.join('%10i'       % v    for v in range(2,len(i)+1)) + '\n')
    of.write('# SED(Lo/A)'  + ' '.join('%s'         % s[v] for v in range(1,len(i)))   + '\n')
    of.write('# Lambda_A '  + ' '.join('%10.4E'     % v    for v in b[1:])             + '\n')
    # Write model sed's to file
    for l in range(len(w)):
        s = f'{w[l]:10.4E}'
        for j in range(1,len(i)):
            s += f'{f[l][int(i[j])]:11.4E}'
        of.write(s+'\n')
    of.close()
    print()

def fdoc(w,t):
    # Lists time scale in fits file
    print()
    print('In file ' + lnam(ifile) + ' there are ' + str('{:d}'.format(len(t))) + ' galaxy SED''s ranging from ' + str('{:.1E}'.format(t[0])) +
          ' to ' + str('{:.1E}'.format(t[len(t)-1])) + ' yr.')
    print('Each SED covers from ' + str('{:.1f}'.format(w[0])) + ' to ' + str('{:.1E}'.format(w[len(w)-1])) + 
           ' A in ' + str('{:d}'.format(len(w))) + ' steps with variable wavelength sampling.')
    print()
    print('The list of time steps follows:')
    print()
    l = len(t)
    n1 = l - (l//7)*7 + 1 ; n2 = l//7 + 1 ; n3 = (l//7)*7 + 1 ; n4 = n3 - 2
    x = ''
    for j in range(8):
        x = x + f'  N    age(yr)      '
    print(x)
    for j in range(1,n1):
        x = ''
        for v in range(j-1,j+n3,n2-1):
            x = x +f'{v+1:3d}: {t[v]:10.3E}     '
        print(x)
    for j in range(n1,n2):
        x = ''
        for v in range(j-1,j+n4,n2-1):
            x = x +f'{v+1:3d}: {t[v]:10.3E}     '
        print(x)
    print()

def bcfits2ised():

    # Writes fortran compatible *.ised2 binary file

    # Open ised2 binary file
    n = ifile.replace('fits','ised2')
    o = FortranFile(n, 'w')
    print('\n               Writing file: ' + n + '    (will take some time...)\n')

    # Define integer parameters
    nb    = len(t)		# number of time steps
    nw    = len(w)		# number of wavelength points
    iseg  = s[0][1]		# number of IMF segments
    jvaz  = s[0][2]		# 
    jseg  = s[0][3]		#
    io    = s[1][0]		# SFR code
    nskip = s[1][1]		# number of dark SEDs to skip
    kdeff = s[1][2]		# BC/CB model code
    sissa = s[1][3]		#  1 if SISSA 2020 stellar tracks used
    nskip = min(0,nskip)	# Check for cb2022 models which contain dark time steps

    # Build integer array
    n = np.array([ -nb, nw, iseg, nskip, sissa, kdeff, io, jvaz, jseg ])
    # Write integer array: 			write (1) nb, nw, iseg, nskip, sissa, kdeff, io, jvaz, jseg
    o.write_record(np.array(n, dtype=np.int32))

    # Write time scale and IMF data:		write (1) tb, ml, mu, (xx(i), lm(i), um(i), baux(i), cn(i), cc(i), i=1, iseg), totm, totn, avs, tau
    ml   = s[1][4]
    mu   = s[1][5]
    totm = s[1][6]
    totn = s[1][7]
    avs  = s[1][8]
    tau  = s[1][9]
    ma = np.array([ ml, mu ])
    mb = np.array([ totm, totn, avs, tau ])
    ia = np.array([ s[0][4], s[0][5], s[0][6], s[0][7], s[0][8], s[0][9] ])
    for j in range (1,iseg):
        ia = np.concatenate((ia, np.array([ s[j][4], s[j][5], s[j][6], s[j][7], s[j][8], s[j][9] ]) ))
    aa = np.concatenate(( t, ma, ia, mb ))
    # Write array aa to file
    o.write_record(np.array(aa, dtype=np.float32))

    # Write wavelength and flux array (f[0] in fits file = wavelength array)
    #                               		write (1) (w(i),   i=1, nw)
    #                               		write (1) (f(i,j), i=1, nw)     for j=1 to nb
    v = [0]*nw
    for j in range(nb+1):
        for i in range(nw):
            v[i] = f[i][j]
        o.write_record(np.array(v, dtype=np.float32))

    # Write physical properties after the spectra
    h = [ 2, 20, 24, 32, 10, 15, 16, 17, 18, 21, 34, 33, 22, 23, 36, 37, 11, 12, 13, 14]
    v = [0]*nb
    for i in range(len(h)):
        for j in range(nb):
            if j < nskip:
                v[j] = 0.
            else:
                if i == 0:
                    v[j] = 10.**(-0.4*(m[j-nskip][h[i]-1]-4.75))
                else:
                    v[j] = p[j-nskip][h[i]-1]
        o.write_record(np.array(v, dtype=np.float32))

    # Close file
    o.close()

def bcfits2ascii():

    # Writes fits file in ASCII format compatible with BC03 release

    # Open ised_ASCII file
    n = ifile.replace('fits','ised_ASCII')
    o = open(n,'w')
    print('\n               Writing file: ' + n + '    (will take some time...)\n')

    # Define integer parameters
    nb    = len(t)		# number of time steps
    nw    = len(w)		# number of wavelength points
    io    = s[1][0]		# SFR code
    nskip = s[1][1]		# number of dark SEDs to skip
    kdeff = s[1][2]		# BC/CB model code
    sissa = s[1][3]		#  1 if SISSA 2020 stellar tracks used
    nskip = min(0,nskip)	# Check for cb2022 models which contain dark time steps

    # Write header
    myheader(o,e,-1)

    # Write time scale
    o.write(str(nb) + ' ' + ' '.join(str(val) for val in t) + '\n')

    # Write wavelength and flux array (f[0] in fits file = wavelength array)
    #                               		write (1) (w(i),   i=1, nw)
    #                               		write (1) (f(i,j), i=1, nw)     for j=1 to nb
    v = [0]*nw
    for j in range(nb+1):
        for i in range(nw):
            v[i] = f[i][j]
        o.write(str(nw) + ' ' + ' '.join(str(val) for val in v) + '\n')

    # Write physical properties after the spectra
    h = [ 2, 20, 24, 32, 10, 15, 16, 17, 18, 21, 34, 33, 22, 23, 36, 37, 11, 12, 13, 14]
    v = [0]*nb
    for i in range(len(h)):
        for j in range(nb):
            if j < nskip:
                v[j] = 0.
            else:
                if i == 0:
                    v[j] = 10.**(-0.4*(m[j-nskip][h[i]-1]-4.75))
                else:
                    v[j] = p[j-nskip][h[i]-1]
        o.write(str(nb) + ' ' + ' '.join(str(val) for val in v) + '\n')

    # Close file
    o.close()

def action(n,k):
    # Writes to file required output
    if k == 0:
        xtab()
    elif k == -2:
        bcfits2ised()
    elif k == -3:
        bcfits2ascii()
    elif (k > 0):
        nsed(n)

def read_bcfits(ifile,k):
    # Reads fits table with BC models
    import warnings
    from   astropy.table import Table
    from   astropy.io import fits
    from   astropy.utils.exceptions import AstropyWarning
    
    # Read file in fits table format
    with fits.open(ifile) as hdul:
        warnings.simplefilter('ignore', category=AstropyWarning)
        # Use Table to read the fits table
        print('Reading file... (will take a few seconds)')
        f = Table.read(hdul,hdu=1)	# hdu = 1 => SED's (luminosity vs. wavelength)
        p = Table.read(hdul,hdu=2)	# hdu = 2 => galaxy physical properties
        m = Table.read(hdul,hdu=3)	# hdu = 3 => photometric magnitude in different bands
        d = Table.read(hdul,hdu=4)	# hdu = 4 => line spectral indices
        t = Table.read(hdul,hdu=5)	# hdu = 5 => time scale for spectral evolution (221 steps)
        s = Table.read(hdul,hdu=6)	# hdu = 6 => miscellaneous data concerning IMF and SFR
        t = t['age-yr']			# time scale for spectral evolution (221 steps)
        w = f['Wavelength']		# wavelength array
        a = 10.**m['logage']		# time scale without pre MS evolution (agrees with t for BC03 and C&B models)
        if a[0] < 10.:
            a[0]=0.			# change first time step to a = 0. It is written as 1 yr in the color files.
        e = hdul[0].header		# ascii table header
        h = hdul[1].header		# SED column header
        for i in range(len(t)+1):
            h[i] = h['TTYPE' + str(i+1)]
    if k==0:
        fdoc(w,t)
    return w,f,t,h,m,d,p,a,s,e

def amup(filew):
    # Checks for MU in file name
    # Mup
    if '_MU600_' in filew:
        m = '/Mu600/'
    elif '_MU300_' in filew:
        m = '/Mu300/'
    elif '_MU010_' in filew:
        m = '/Mu010/'
    else:
        m = '/Mu100/'
    return m

def sspdir(filew):
    # Build default directory for this file name
    # Model:
    q = True
    d = False
    if 'bc2003' in filew:
        mdir = 'BC03'
    elif 'cb2003' in filew:
        mdir = 'CB03'
    elif 'cb2007' in filew:
        mdir = 'CB07'
    elif 'bc2019' in filew:
        mdir = 'BC19'
    elif 'cb2019' in filew:
        mdir = 'CB19'
        d = True
    elif 'bc2022' in filew:
        mdir = 'BC22'
    elif 'cb2022' in filew:
        mdir = 'CB22'
    else:
        q = False

    # IMF:
    if '_chab_' in filew:
        mdir = mdir + '_chabrier'
    elif '_kroup_' in filew:
        mdir = mdir + '_kroupa'
    elif '_salp_' in filew:
        mdir = mdir + '_salpeter'
    elif '_v0p30_' in filew:
        mdir = mdir + '_vazdekis_0.30'
    elif '_v0p80_' in filew:
        mdir = mdir + '_vazdekis_0.80'
    elif '_v1p00_' in filew:
        mdir = mdir + '_vazdekis_1.00'
    elif '_v1p30_' in filew:
        mdir = mdir + '_vazdekis_1.30'
    elif '_v1p50_' in filew:
        mdir = mdir + '_vazdekis_1.50'
    elif '_v1p80_' in filew:
        mdir = mdir + '_vazdekis_1.80'
    elif '_v2p00_' in filew:
        mdir = mdir + '_vazdekis_2.00'
    elif '_v2p30_' in filew:
        mdir = mdir + '_vazdekis_2.30'
    elif '_v2p80_' in filew:
        mdir = mdir + '_vazdekis_2.80'
    elif '_v3p30_' in filew:
        mdir = mdir + '_vazdekis_3.30'
    else:
        q = False

    # Mup
    if d:
        m = amup(filew)
    else:
        m = '/'

    if q:
        mdir = os.environ.get('glxssp') + '/' + mdir + m
    else:
        mdir = ''
    return mdir

def fcheck(file):
    # Checks if file exists in various directories
    file = lnam(file)
    if str(file).find('.fits') < 0:
        file = file + '.fits'
    glxext = sspdir(file)
    while True:
        if os.path.isfile(os.getcwd() + '/' + file):
            file = os.getcwd() + '/' + file
           #print('cwd','   ',file)
            break
        elif os.path.isfile(os.environ.get('glxssp') + '/' + file):
            file = os.environ.get('glxssp') + '/' + file
           #print('ssp','   ',file)
            break
        elif os.path.isfile(os.environ.get('glxout') + '/' + file):
            file = os.environ.get('glxout') + '/' + file
           #print('out','   ',file)
            break
        elif os.path.isfile(glxext + file):
            file = glxext + file
           #print('ext','   ',file)
            break
        else:
           #file = input('File ' + rp(file) + ' does not exist. Enter new file name = ')
            file = input('File ' +    file  + ' does not exist. Enter new file name = ')
    return file

def settings():
    # Read environment variables
    from dotenv import load_dotenv
   #load_dotenv(os.getenv('GALAXEV') + '/pylib/')
    load_dotenv()

# Execute selected option - Command Line version

# Read command line
g = sys.argv
k = len(g)

settings()

if k < 2:
    print('Usage:  bcfits2txt.py model_file.fits')
else:
    # Read fits file
    ifile = g[1]
    ifile = fcheck(ifile)
    w,f,t,h,m,d,p,a,s,e = read_bcfits(ifile,0)

    if k == 2:
        print ('Select records to extract:')
        print ('                            t1 t2 t3... (age in Gyr of SED''s)')
        print ('                           -N1 N2 N3... (record number of SED''s. Note - sign in front of N1)')
        print ('                            N1:N2       (all SED''s with N1 <= N <= N2)')
        print ('                            a           (all SED''s in fits file)')
        print ('                            t           (write tables with different properties of SSP model)')
        print ('                            i           (create .ised fortran readable binary file compatible with BC03 model release)')
        print ('                            c           (create .ised_ASCII text file compatible with BC03 model release)')
        print ('                            h           (enter your choice in the command line, examples)')
        print ('                            ENTER       (exits loop)')
        while True:
            # Ask for record to extract
            n,k = getc('                    Choice: ')
            if len(n) <= 0:
                break
            action(n,k)
    else:
        # command line mode
        n,k = getb(g[2])
        action(n,k)

import os
import math
import builtins as bt
import basics.basics as bs

def fillin(s,u):
    for i in range(u-len(s)-1):
        s = s + ' '
    s = s + 'I'
    return s

def sedid(a,n,u,f):

    def oprint(s,f):
        if f=='s':
            print(s)
        else:
            f.write(s+'\n')

    if n == 2022:
        # BaSeL 3.1 stellar library model
        s = 'I      S.E.D.: UV to IR -  91 A to 360E6 A: BaSeL, Rauch and Aringer et al. stellar models + IRTF stellar library +'    ; s = fillin(s,u) ; oprint(a+s,f)
        s = 'I                                           Dusty code models for TP-AGB stars'					     ; s = fillin(s,u) ; oprint(a+s,f)
    elif n == 7124:
        # Stelib stellar library model
        s = 'I      S.E.D.: Visible - 3322 A to  8750 A: Stelib stellar library'						     ; s = fillin(s,u) ; oprint(a+s,f)
        s = 'I                   UV -   91 A to  3310 A: BaSeL  stellar library'       						     ; s = fillin(s,u) ; oprint(a+s,f)
        s = 'I                   IR - 8770 A to 360E6 A: BaSeL and Aringer et al. stellar models + IRTF stellar library +'	     ; s = fillin(s,u) ; oprint(a+s,f)
        s = 'I                                           Dusty code models for TP-AGB stars'					     ; s = fillin(s,u) ; oprint(a+s,f)
    elif n == 16902:
        # Miles stellar library model
        s = 'I      S.E.D.: Visible - 3540 A to  7350 A: Miles stellar library'                                                	     ; s = fillin(s,u) ; oprint(a+s,f)
        s = 'I                        7350 A to  9400 A: IndoUS stellar library'                                               	     ; s = fillin(s,u) ; oprint(a+s,f)
        s = 'I                   UV -    5 A to  3540 A: Tlusty, Martins et al, UVBlue, Rauch, WMBasic and PoWR stellar models'	     ; s = fillin(s,u) ; oprint(a+s,f)
        s = 'I                                           resampled every 0.9 A using the flux preseving procedure by Carnall (2017)' ; s = fillin(s,u) ; oprint(a+s,f)
        s = 'I                   IR - 9410 A to 360E6 A: BaSeL and Aringer et al. stellar models + IRTF stellar library +'     	     ; s = fillin(s,u) ; oprint(a+s,f)
        s = 'I                                           Dusty code models for TP-AGB stars'					     ; s = fillin(s,u) ; oprint(a+s,f)
    elif n == 46226:
        # Miles stellar library model
        s = 'I      S.E.D.: Visible - 3540 A to  7350 A: Miles stellar library'                                                	     ; s = fillin(s,u) ; oprint(a+s,f)
        s = 'I                        7350 A to  9400 A: IndoUS stellar library'                                               	     ; s = fillin(s,u) ; oprint(a+s,f)
        s = 'I                   UV -    5 A to  3540 A: Tlusty, Martins et al, UVBlue, Rauch, WMBasic and PoWR stellar models'	     ; s = fillin(s,u) ; oprint(a+s,f)
        s = 'I                                           resampled every 0.1 A using the flux preseving procedure by Carnall (2017)' ; s = fillin(s,u) ; oprint(a+s,f)
        s = 'I                   IR - 9410 A to 360E6 A: BaSeL and Aringer et al. stellar models + IRTF stellar library +'     	     ; s = fillin(s,u) ; oprint(a+s,f)
        s = 'I                                           Dusty code models for TP-AGB stars'					     ; s = fillin(s,u) ; oprint(a+s,f)
    elif n == 24766:
        # IndoUS stellar library model
        s = 'I      S.E.D.: Visible - 3540 A to  9400 A: IndoUS stellar library'                                               	     ; s = fillin(s,u) ; oprint(a+s,f)
        s = 'I                   UV -    5 A to  3540 A: Tlusty, Martins et al, UVBlue, Rauch, WMBasic and PoWR stellar models'	     ; s = fillin(s,u) ; oprint(a+s,f)
        s = 'I                   IR - 9410 A to 360E6 A: BaSeL and Aringer et al. stellar models + IRTF stellar library +'     	     ; s = fillin(s,u) ; oprint(a+s,f)
        s = 'I                                           Dusty code models for TP-AGB stars'                                   	     ; s = fillin(s,u) ; oprint(a+s,f)
    s = 'I'
    s = fillin(s,u)
    oprint(a+s,f)

def idheader(e,n):
    # Shows model id (header) on the screen
    a = '     '
    i = 0
    print()
    for k in range(1000):
        if 'END' in e['HEADER' + str(i+1)]:
           break
        else:
            s = e['HEADER' + str(i+1)] + e['HEADER' + str(i+2)]
            s = s.replace('|',' ')
            if k == 0:
                d = 'I      BUILT:  ' + s[len(s)-24:len(s)] 		# date
            if k==1:
                u = len(s)
            if k==3:
                j = s.find('=') + 2
                l = s[j:].find(' ')
                s = s[:j-1] + ' ' + bs.bp(s[j:j+l]) + s[j+l+1:len(s)]
                s = fillin(s,u)
                s = s.replace('II',' I')
           #if k==3 and len(s) > u:
           #    s = s.replace('  I','I')
            if k==5:
                d = fillin(d,u)
                print(a+d)
                d = 'I'
                d = fillin(d,u)
                print(a+d)
            if k==7:
                sedid(a,n,u,'s')
            if k>0:
                print(a+s)
            if 'totals' in s:
                print(a + 'I                                                                                                                        I')
                print(a + 'I                              (C) 1995-2023 G. Bruzual A. & S. Charlot - All Rights Reserved                            I')
                print(a + 'I------------------------------------------------------------------------------------------------------------------------I')
                break
            i=i+2

def ytab(ifile,w,t,e,m,d,p,v,a):
    # Writes header in BC tables
    global zu

    def dark(of):
        # Writes to file time steps in which the galaxy is still dark
        for i in range(len(t)-len(a)):
            if i==0:
                s = '#%9.6f' % 0.00 + '    Dark galaxy'
            else:
                s = '#%9.6f' % math.log10(t[i]) + '    Dark galaxy'
            of.write(s + '\n')

    x = ['.0VEGAmag', '.1VEGAmag', '.2ABmag', '.3ABmag', '.4lsindx', '.5lsindx', '.6lsindx', '.7physp', '.8physp', '.9physp', '.fwage_rf']
    f = bs.lnam(ifile.replace('.fits',''))
    o = (os.environ.get('glxout') + '/' + f).replace('//','/')
    j = ' '
    print(27*j,'Creating files: ' + o + '.0VEGAmag ...')
    print(50*j,'+ (.1VEGAmag, .2ABmag, .3ABmag, .4lsindx, .5lsindx, .6lsindx, .7physp, .8physp, .9physp, .fwage_rf)')
    print()
    a = m['logage']				# time scale without pre MS evolution (agrees with t for BC03 and C&B models)
    if len(v[0]) >= 3:
        # Output from csp_galaxev, includes .fwage_rf
        lx = len(x)
        zu = v[0][len(v[0])-1]
    else:
        # Fits file does not include .fwage_rf
        lx = len(x) - 1
    for l in range(lx):
        of = o + x[l]
        of = open(of,'w')
        myheader2(of,w,t,e,l)

        if l == 0:
            # Vega mags
            h = [ 2, 14, 23, 22, 32, 43, 49, 44, 53, 59, 69, 71, 61, 63, 70, 66]
            u = [0]*len(h)
            dark(of)
            for r in range(len(a)):
                for j in range(len(h)):
                    u[j] = m[r][h[j]-1]
                s = '%10.6f' % a[r] + ' ' + ' '.join('%9.4f' % v for v in u)
                of.write(s + '\n')
            of.close()

        if l == 1:
            # Vega mags
            h = [ 2, 32, 58, 65, 68, 72, 73, 74, 75, 76, 78, 79, 81, 77, 80, 82]
            u = [0]*len(h)
            dark(of)
            for r in range(len(a)):
                for j in range(len(h)):
                    u[j] = m[r][h[j]-1]
                s = '%10.6f' % a[r] + ' ' + ' '.join('%9.4f' % v for v in u)
                of.write(s + '\n')
            of.close()

        if l == 2:
            # AB mags
            h = [13, 24, 35, 45, 55, 16, 15, 27, 26, 36, 39, 46, 47, 54, 56, 64, 67, 3, 5, 83, 84, 85]
            u = [0]*len(h)
            dark(of)
            for r in range(len(a)):
                for j in range(len(h)):
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
            u = [0]*len(h)
            dark(of)
            for r in range(len(a)):
                for j in range(len(h)):
                    u[j] = m[r][h[j]-1]
                s = '%10.6f' % a[r] + ' ' + ' '.join('%9.4f' % v for v in u)
                of.write(s + '\n')
            of.close()

        if l == 4:
            # Lick indices
            u = [0]*len(h)
            dark(of)
            for r in range(len(a)):
                for j in range(1,22):
                    u[j-1] = d[r][j]
                s = '%10.6f' % a[r] + ' ' + ' '.join('%8.4f' % v for v in u)
                of.write(s + '\n')
            of.close()

        if l == 5:
            # Other indices
            u = [0]*len(h)
            dark(of)
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
            u = [0]*len(h)
            dark(of)
            for r in range(len(a)):
                for j in range(36,57):
                    u[j-36] = d[r][j]
                s = '%10.6f' % a[r] + ' ' + ' '.join('%10.4f' % v for v in u)
                of.write(s + '\n')
            of.close()

        if l == 7:
            # Physical properties
            h = [ 9, 7, 8, 6, 2, 3, 4, 5, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
            u = [0]*len(h)
            dark(of)
            for r in range(len(a)):
                for j in range(len(h)):
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
            u = [0]*len(h)
            dark(of)
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
            u = [0]*len(h)
            dark(of)
            for r in range(len(a)):
                for j in range(1):
                    u[j] = m[r][h[j]-1]
                for j in range(1,7):
                    u[j] = p[r][h[j]-1]
                s = '%10.6f' % a[r] + '%10.4f' % u[0] + ' ' + ' '.join('%12.4E' % z for z in u[1:8])
                of.write(s + '\n')
            of.close()

        if l == 10:
            # Flux weighted age in galaxy restframe at z = zu
            h = [ 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
            u = [0]*len(h)
            dark(of)
            for r in range(len(a)):
                for j in range(len(h)):
                    u[j] = v[r][h[j]]
                s = '%10.6f' % a[r] + ' ' + ' '.join('%12.4E' % c for c in u[0:7]) + ' ' + ' '.join('%11.5f' % c for c in u[7:len(h)])
                of.write(s + '\n')
            of.close()

def myheader2(f,w,t,e,l):
    # Writes header in BC tables
    a = '#          '
    h = ['3', '3', '3', '3', '5', '5', '5', '3', '3', '3', '3']
    i = 0
    for k in range(1000):
        if 'END' in e['HEADER' + str(i+1)]:
           break
        else:
            s = e['HEADER' + str(i+1)] + e['HEADER' + str(i+2)]
            s = s.replace('|',' ')
            if k==1:
                u = len(s)
            if k==7:
                sedid(a,len(w),u,f)
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
        h4 = '# log-age    HdeltaA  HgammaA  HdeltaF  HgammaF    B4_GC    B4_VN   CaII8498   CaII8542   CaII8662   MgI8807   H8_3889   H9_3835  H10_3798     BH-HK'
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

    elif l == 10:
        # Header for file .fwage_rf
        h1 = '#   (1)         (2)          (3)          (4)          (5)          (6)          (7)          (8)           (9)        (10)        (11)        (12)        (13)        (14)        (15)'
        h2 = "#            <-------------------------------------------------------- computed after red shifting sed's by z = " + f'{zu:6.3f} ----------------------------------------------------------------->'
        h3 = '# log_t(yr)    MWA(yr)     uFWA(yr)     gFWA(yr)     rFWA(yr)     iFWA(yr)     zFWA(yr)     KFWA(yr)     MWLA(yr)   uFWLA(yr)   gFWLA(yr)   rFWLA(yr)   iFWLA(yr)   zFWLA(yr)   KFWLA(yr)'

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


def pid(f,w,t):
    # prints model basic data
    print()
    print(" In file " + bs.bp(bs.lnam(f)) + " there are " + str('{:d}'.format(len(t))) + " galaxy SED's, ranging from " + str('{:.1E}'.format(t[0])) +
        "  to " + str('{:.1E}'.format(t[len(t)-1])) + " yr")
    print(" Each SED covers the wavelength range from " + str('{:.1f}'.format(w[0])) + " to " + str('{:.1E}'.format(w[len(w)-1])) +
        "  A in " + str('{:d}'.format(len(w))) + " steps")
    print()

def fdoc(t):
    # Lists time scale in fits file
    print(' The list of time steps follows:')
    print()
    l = len(t)
    n1 = l - (l//7)*7 + 1 ; n2 = l//7 + 1 ; n3 = (l//7)*7 + 1 ; n4 = n3 - 2
    x = ''
    for j in range(8):
       #x = x + f'  N     age(yr)        '
        x = x + f'  N     age(yr)      '
    print(x)
    for j in range(1,n1):
        x = ''
        for v in range(j-1,j+n3,n2-1):
           #x = x +f'{v+1:3d}: {t[v]:13.6E}     '
            x = x +f'{v+1:3d}: {t[v]:12.5E}    '
        print(x)
    for j in range(n1,n2):
        x = ''
        for v in range(j-1,j+n4,n2-1):
           #x = x +f'{v+1:3d}: {t[v]:13.6E}     '
            x = x +f'{v+1:3d}: {t[v]:12.5E}    '
        print(x)
    print()

def multhead(k,f,t,w,e):
    # Prints various id headers
    if k == 1:
        idheader(e,len(w))
        pid(f,w,t)
        fdoc(t)
    elif k == 2:
       #idheader(e,len(w))
       #pid(f,w,t)
        fdoc(t)
   #else:
   #    pid(f,w,t)

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

def printmag(t,z,dm,tl,tz,av,ab,fj,o,i):
    if i == 0:
        a = t*1.E9
        s1 = '%10.3E' % a + ' ' + '%9.4f' % z + ' ' + '%9.4f' % dm + ' ' + '%7.2f' % av + ' ' + ' '.join('%9.4f' % c for c in ab)
        if not o:
            print(s1)
        else:
            s2 = '%10.3E' % a + ' ' + '%9.4f' % z + ' ' + '%9.4f' % dm + ' ' + '%7.2f' % av + ' ' + ' '.join('%12.4E' % c for c in fj)
            o[1].write(s1 + '\n')
            o[2].write(s2 + '\n')
    else:
        s1 = f'{z:8.4f}{tl:10.5f}{tz:10.5f}{dm:9.4f}{av:7.2f}' + ' ' + ' '.join('%9.4f' % c for c in ab)
        o[1].write(s1 + '\n')

def hdusr(i):
    # Builds header for user selected filters
    l = []
    for j in range(len(i)):
        f = 'FLT' + f'{i[j]:3d}'
        f = f.replace(' ','0')
        l.append(f)
    return l

def opfiles(inp,z,av):
    # Open output files and writes header. Used by program rf_phot
    fds  = inp[1]
    if fds=='john':
        fds='johnson'
    if '-o' in inp[len(inp)-1]:
        file = inp[0]
        n = file.replace('fits','rf_phot.' + fds)
        p = True
        n = n + ('.z' + f'{z:6.3f}' + '.Av' + f'{av:6.3f}').replace(' ','')
    elif '-o' in inp[len(inp)-2]:
        n = inp[len(inp)-1]
        p = True
    else:
        p = False
    if p:
        if fds=='jplus':
            f = os.environ.get('glxout') + '/' + n + '.ABmag'
            g = f.replace('ABmag','Fjansky')
            o = open(f,'w')
            q = open(g,'w')
            o.write('#   t_yr        z        dm        Av      uJAVA     J0378     J0395     J0410     J0430     gSDSS     J0515     rSDSS     J0660     iSDSS     J0861     zSDSS\n')
            q.write('#   t_yr        z        dm        Av      uJAVA        J0378        J0395        J0410        J0430        gSDSS        J0515        rSDSS        J0660        iSDSS        J0861        zSDSS\n')
        elif fds=='jpas':
            v=' '
            f = os.environ.get('glxout') + '/' + n + '.ABmag'
            g = f.replace('ABmag','Fjansky')
            o = open(f,'w')
            q = open(g,'w')
            o.write('#' + 31*v + 'Column:' + '   '.join('%7d' % c for c in range(4,64)) + '\n')
            o.write('#' + 31*v + 'Filter:' + '   '.join('%7d' % c for c in range(1,61)) + '\n')
            o.write('#   t_yr        z        dm        Av'
                    +6*v+'uJAVA'+7*v+'u'    +7*v+'J0378'+5*v+'J0390'+5*v+'J0400'+5*v+'J0410'+5*v+'J0420'+5*v+'J0430'+5*v+'J0440'+5*v+'J0450'+5*v+'J0460'+5*v+'J0470' 
                    +5*v+'J0480'+5*v+'gSDSS'+5*v+'J0490'+5*v+'J0500'+5*v+'J0510'+5*v+'J0520'+5*v+'J0530'+5*v+'J0540'+5*v+'J0550'+5*v+'J0560'+5*v+'J0570'+5*v+'J0580'
                    +5*v+'J0590'+5*v+'J0600'+5*v+'J0610'+5*v+'J0620'+5*v+'rSDSS'+5*v+'J0630'+5*v+'J0640'+5*v+'J0650'+5*v+'J0660'+5*v+'J0670'+5*v+'J0680'+5*v+'J0690'
                    +5*v+'J0700'+5*v+'J0710'+5*v+'J0720'+5*v+'J0730'+5*v+'J0740'+5*v+'J0750'+5*v+'J0760'+5*v+'iSDSS'+5*v+'J0770'+5*v+'J0780'+5*v+'J0790'+5*v+'J0800'
                    +5*v+'J0810'+5*v+'J0820'+5*v+'J0830'+5*v+'J0840'+5*v+'J0850'+5*v+'J0860'+5*v+'J0870'+5*v+'J0880'+5*v+'J0890'+5*v+'J0900'+5*v+'J0910'+5*v+'J1007\n')
            q.write('#' + 31*v + 'Column:      4  ' + '  '.join('%11d' % c for c in range(5,64)) + '\n')
            q.write('#' + 31*v + 'Filter:      1  ' + '  '.join('%11d' % c for c in range(2,61)) + '\n')
            q.write('#   t_yr        z        dm        Av'
                    +6*v+'uJAVA'+9*v+' u  ' +8*v+'J0378'+8*v+'J0390'+8*v+'J0400'+8*v+'J0410'+8*v+'J0420'+8*v+'J0430'+8*v+'J0440'+8*v+'J0450'+8*v+'J0460'+8*v+'J0470' 
                    +8*v+'J0480'+8*v+'gSDSS'+8*v+'J0490'+8*v+'J0500'+8*v+'J0510'+8*v+'J0520'+8*v+'J0530'+8*v+'J0540'+8*v+'J0550'+8*v+'J0560'+8*v+'J0570'+8*v+'J0580'
                    +8*v+'J0590'+8*v+'J0600'+8*v+'J0610'+8*v+'J0620'+8*v+'rSDSS'+8*v+'J0630'+8*v+'J0640'+8*v+'J0650'+8*v+'J0660'+8*v+'J0670'+8*v+'J0680'+8*v+'J0690'
                    +8*v+'J0700'+8*v+'J0710'+8*v+'J0720'+8*v+'J0730'+8*v+'J0740'+8*v+'J0750'+8*v+'J0760'+8*v+'iSDSS'+8*v+'J0770'+8*v+'J0780'+8*v+'J0790'+8*v+'J0800'
                    +8*v+'J0810'+8*v+'J0820'+8*v+'J0830'+8*v+'J0840'+8*v+'J0850'+8*v+'J0860'+8*v+'J0870'+8*v+'J0880'+8*v+'J0890'+8*v+'J0900'+8*v+'J0910'+8*v+'J1007\n')
        elif fds== 'sdss':
            f = os.environ.get('glxout') + '/' + n + '.ABmag'
            g = f.replace('ABmag','Fjansky')
            o = open(f,'w')
            q = open(g,'w')
            o.write('#   t_yr        z        dm        Av       u         g         r         i         zz\n')
            q.write('#   t_yr        z        dm        Av        u            g            r            i            z\n')
        elif fds== 'johnson':
            f = os.environ.get('glxout') + '/' + n + '.VEGAmag'
            g = f.replace('VEGAmag','Nphotons')
            o = open(f,'w')
            q = open(g,'w')
            o.write('#   t_yr        z        dm        Av       U        B2        B3         V        Rc        Ic         R         I         J         K         L        PalJ      PalH      PalK       K\'\n')
            q.write('#   t_yr        z        dm        Av        U           B2           B3            V           Rc           Ic            R            I            J            K            L           PalJ         PalH         PalK          K\'\n')
        elif fds== 'stmag':
            f = os.environ.get('glxout') + '/' + n + '.STmag'
            g = f.replace('STmag','STflux')
            o = open(f,'w')
            q = open(g,'w')
            o.write('#                                        <----------------------------------------------------------- HST ACS STmag ------------------------------------------------------------->   <--------------------------- JWST NIRCAM STmag ----------------------------->   <---------------------------------- JWST MIRI STmag ---------------------------------->\n')
            o.write('#   t_yr        z        dm        Av     F150LP    F165LP     F220w     F250w     F336W     F410w     F438W     F475w     F547M     F555w     F606w     F625w     F775w     F814w     F070W     F090W     F115W     F150W     F200W     F277W     F356W     F444W    F0560W    F0770W    F1000W    F1130W    F1280W    F1500W    F1800W    F2100W    F2550W\n')
            q.write('#                                        <-------------------------------------------------------------------------------- HST ACS STflux --------------------------------------------------------------------------------->   <--------------------------------------- JWST NIRCAM STflux ---------------------------------------->   <------------------------------------------------ JWST MIRI STflux ---------------------------------------------->\n')
            q.write('#   t_yr        z        dm        Av      F150LP       F165LP       F220w        F250w        F336W        F410w        F438W        F475w        F547M        F555w        F606w        F625w        F775w        F814w        F070W        F090W        F115W        F150W        F200W        F277W        F356W        F444W        F0560W       F0770W       F1000W       F1130W       F1280W       F1500W       F1800W       F2100W       F2550W\n')
        elif 'user' in fds:
            lu =hdusr(bt.usr)
            f = os.environ.get('glxout') + '/' + n
            if 'ABmag' in f:
                f = f.replace('.ABmag','') + '.ABmag'
                g = f.replace('ABmag','Fjansky')
            elif 'STmag' in f:
                f = f.replace('.STmag','') + '.STmag'
                g = f.replace('STmag','STflux')
            elif 'VEGAmag' in f:
                f = f.replace('.VEGAmag','') + '.VEGAmag'
                g = f.replace('VEGAmag','Nphotons')
            o = open(f,'w')
            q = open(g,'w')
            o.write('#   t_yr        z        dm        Av ' + ''.join('%10s' % c for c in (lu)) + '\n')
            q.write('#   t_yr        z        dm       Av'   + ''.join('%13s' % c for c in (lu)) + '\n')
        w = True
        return w,o,q,f,g
    else:
        if fds == 'jplus':
            print('#   t_yr        z        dm        Av      uJAVA     J0378     J0395     J0410     J0430     gSDSS     J0515     rSDSS     J0660     iSDSS     J0861     zSDSS')
        elif fds == 'jpas':
            print('#   t_yr        z        dm        Av      uJAVA       u       J0378     J0390     J0400     J0410     J0420     J0430     J0440     J0450     J0460     J0470'+ 
                                                        '     J0480     gSDSS     J0490     J0500     J0510     J0520     J0530     J0540     J0550     J0560     J0570     J0580'+
                                                        '     J0590     J0600     J0610     J0620     rSDSS     J0630     J0640     J0650     J0660     J0670     J0680     J0690'+
                                                        '     J0700     J0710     J0720     J0730     J0740     J0750     J0760     iSDSS     J0770     J0780     J0790     J0800'+
                                                        '     J0810     J0820     J0830     J0840     J0850     J0860     J0870     J0880     J0890     J0900     J0910     J1007')
        elif fds == 'sdss':
            print('#   t_yr        z        dm        Av       u         g         r         i         z')
        elif fds == 'johnson':
            print('#   t_yr        z        dm        Av       U        B2        B3         V        Rc        Ic         R         I         J         K         L        PalJ      PalH      PalK       K\'')
        elif fds == 'stmag':
            print('#                                        <----------------------------------------------------------- HST ACS STmag ------------------------------------------------------------->   <--------------------------- JWST NIRCAM STmag ----------------------------->   <---------------------------------- JWST MIRI STmag ---------------------------------->')
            print('#   t_yr        z        dm        Av     F150LP    F165LP     F220w     F250w     F336W     F410w     F438W     F475w     F547M     F555w     F606w     F625w     F775w     F814w     F070W     F090W     F115W     F150W     F200W     F277W     F356W     F444W    F0560W    F0770W    F1000W    F1130W    F1280W    F1500W    F1800W    F2100W    F2550W')
        elif 'user' in fds:
            print('#   t_yr        z        dm        Av ' + ''.join('%10s' % c for c in (lu)))
        w = False
        return w

def oqfiles(inp,av):
    # Open output files and writes header. Used by program of_phot
    fds  = inp[1]
    if fds=='john':
        fds='johnson'

    # Cosmological parameters
    h, q, clambda, tu, ttg, zf, omega, omega_lambda = bt.cosmol2
    l='# Cosmology: Ho = '+f'{h:4.0f}'+', Omega ='+f'{omega:5.2f}'+', Omega_lambda ='+f'{omega_lambda:5.2f}'+', qo = '+f'{q:7.4f}'+', Lambda ='+f'{clambda:10.3E}'+', tu = '+f'{tu:5.2f} Gyr'+', tg = '+f'{ttg:5.2f} Gyr'+',  zf = '+f'{zf:6.2f}\n'

    # Open file and write header
    if '-o' in inp[len(inp)-2]:
        n = inp[len(inp)-1]
    else:
        file = inp[0]
        n = file.replace('fits','of_phot.' + fds)
        n = n + ('.Av' + f'{av:6.3f}').replace(' ','')

    # Build header
    v=' '
    if fds=='jplus':
        f = os.environ.get('glxout') + '/' + n + '.ABmag'
        o = open(f,'w')
        o.write(l)
        o.write('#\n')
        o.write('#    z    LTT(Gyr)   tz(Gyr)    dm       Av'
                +6*v+'uJAVA     J0378     J0395     J0410     J0430     gSDSS     J0515     rSDSS     J0660     iSDSS     J0861     zSDSS\n')
    elif fds=='jpas':
        f = os.environ.get('glxout') + '/' + n + '.ABmag'
        o = open(f,'w')
        o.write(l)
        o.write('#\n')
        o.write('#' + 37*v + 'Column:' + '   '.join('%7d' % c for c in range(5,65)) + '\n')
        o.write('#' + 37*v + 'Filter:' + '   '.join('%7d' % c for c in range(1,61)) + '\n')
        o.write('#    z    LTT(Gyr)   tz(Gyr)    dm       Av'
                +6*v+'uJAVA'+7*v+'u'    +7*v+'J0378'+5*v+'J0390'+5*v+'J0400'+5*v+'J0410'+5*v+'J0420'+5*v+'J0430'+5*v+'J0440'+5*v+'J0450'+5*v+'J0460'+5*v+'J0470' 
                +5*v+'J0480'+5*v+'gSDSS'+5*v+'J0490'+5*v+'J0500'+5*v+'J0510'+5*v+'J0520'+5*v+'J0530'+5*v+'J0540'+5*v+'J0550'+5*v+'J0560'+5*v+'J0570'+5*v+'J0580'
                +5*v+'J0590'+5*v+'J0600'+5*v+'J0610'+5*v+'J0620'+5*v+'rSDSS'+5*v+'J0630'+5*v+'J0640'+5*v+'J0650'+5*v+'J0660'+5*v+'J0670'+5*v+'J0680'+5*v+'J0690'
                +5*v+'J0700'+5*v+'J0710'+5*v+'J0720'+5*v+'J0730'+5*v+'J0740'+5*v+'J0750'+5*v+'J0760'+5*v+'iSDSS'+5*v+'J0770'+5*v+'J0780'+5*v+'J0790'+5*v+'J0800'
                +5*v+'J0810'+5*v+'J0820'+5*v+'J0830'+5*v+'J0840'+5*v+'J0850'+5*v+'J0860'+5*v+'J0870'+5*v+'J0880'+5*v+'J0890'+5*v+'J0900'+5*v+'J0910'+5*v+'J1007\n')
    elif fds== 'sdss':
        f = os.environ.get('glxout') + '/' + n + '.ABmag'
        o = open(f,'w')
        o.write(l)
        o.write('#\n')
        o.write('#    z    LTT(Gyr)   tz(Gyr)    dm       Av       u         g         r         i         zz\n')
    elif fds== 'johnson':
        f = os.environ.get('glxout') + '/' + n + '.VEGAmag'
        o = open(f,'w')
        o.write(l)
        o.write('#\n')
        o.write('#    z    LTT(Gyr)   tz(Gyr)    dm       Av       U        B2        B3         V        Rc        Ic         R         I         J         K         L        PalJ      PalH      PalK       K\'\n')
    elif fds== 'stmag':
        f = os.environ.get('glxout') + '/' + n + '.STmag'
        o = open(f,'w')
        o.write(l)
        o.write('#\n')
        o.write('#                                              <----------------------------------------------------------- HST ACS STmag ------------------------------------------------------------->   <--------------------------- JWST NIRCAM STmag ----------------------------->   <---------------------------------- JWST MIRI STmag ---------------------------------->\n')
        o.write('#    z    LTT(Gyr)   tz(Gyr)    dm       Av     F150LP    F165LP     F220w     F250w     F336W     F410w     F438W     F475w     F547M     F555w     F606w     F625w     F775w     F814w     F070W     F090W     F115W     F150W     F200W     F277W     F356W     F444W    F0560W    F0770W    F1000W    F1130W    F1280W    F1500W    F1800W    F2100W    F2550W\n')
    elif 'user' in fds:
        lu =hdusr(bt.usr)
        f = os.environ.get('glxout') + '/' + n
        if 'ABmag' in f:
            f = f.replace('.ABmag','') + '.ABmag'
        elif 'STmag' in f:
            f = f.replace('.STmag','') + '.STmag'
        elif 'VEGAmag' in f:
            f = f.replace('.VEGAmag','') + '.VEGAmag'
        o = open(f,'w')
        o.write('#    z    LTT(Gyr)   tz(Gyr)    dm       Av ' + ''.join('%10s' % c for c in (lu)) + '\n')
    w = True
    return w,o,o,f

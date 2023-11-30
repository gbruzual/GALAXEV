import builtins as bt
import numpy as np
from   astropy.table import Table
import bcfilt.bcfilt as fl
import bcindx.bcindx as ix

def opent(i):
    # Open temporary files needed by BC model codes
    # global td1,td2,td3,td4,td5,td6,td7   #  , o1,o2,o3,o4,o5,o6

    # Open files
   #o1 = open(os.environ.get('glxout') + '/' + f + '.sed','w')
   #o2 = open(os.environ.get('glxout') + '/' + f + '.physical_properties','w')
   #o3 = open(os.environ.get('glxout') + '/' + f + '.photometry','w')
   #o4 = open(os.environ.get('glxout') + '/' + f + '.indices','w')
   #o5 = open(os.environ.get('glxout') + '/' + f + '.time_tsteps','w')
   #o6 = open(os.environ.get('glxout') + '/' + f + '.misc','w')

    # Write headers to files/tables

    # Init table to store model SEDs
    bt.td1 = Table()

    # Physical property file/table
    kp  = [1, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22, 23, 31, 33, 34, 35, 36]
    s   = ['logage_yr', 'NLy', 'NHeI', 'NHeII', 'NHeII/NLy', 'B912', 'B4_VN', 'B4_SDSS', 'B4000', 'SNe/yr', 'PISNe/yr', 'RegSNe/yr', \
           'IaSNe/yr', 'FailedSNe/yr', 'PNBR/yr', 'N(BH)', 'N(NS)', 'N(WD)', 'Lx(Lo)', 'M*_liv', 'M_remnants', 'M_ret_gas', 'M_galaxy', \
           'SFR/yr', 'M*_liv+M_rem', 'M*_tot/Lb', 'M*_tot/Lv', 'M*_tot/Lk', 'M*_liv/Lb', 'M*_liv/Lv', 'M*_liv/Lk', "b(t)*'s/yr", 'B(t)/yr/Lo', \
           'Turnoff_mass', 'BPMS/BMS', 'Loss Rate', 'Prod. Rate' ]
    bt.td2 = Table(names=s, dtype=['f4']*len(s))

    # Photometry file/table
    # Define filters to compute magnitudes (kf > 0, VEGA system; kf < 0, AB system)
    kf  = [-139, -219, -140, -200, -242, -243, -220, -221, -244, -201, -120,   12, -266, -237, -222, -202, -223, -203, -245,   14, \
             13, -121, -224, -267, -238, -247, -204, -225, -246,   15, -226, -248, -122, -239, -249, -227, -268, -250, -251, -252, \
             84,   32, -123, -255, -269, -228,   85, -253, -229, -205,   33, -241, -124, -270, -197,  125,   34, -198,   55, -199, \
             56, -126,  126,   88, -256,  127,   35,   57,   36,  128,  129,  130,  131,   71,  132,   72,   73,  133,   74,  134, -235 ]
    s   = ['logage_yr', 'Mbol', 'FUV_GALEX_AB', 'ACSWFC_F220w_AB', 'NUV_GALEX_AB', 'WFC3_F225W_AB', 'UVIS1_f225w_AB', \
           'UVIS1_f275w_AB', 'ACSWFC_F250w_AB', 'ACSWFC_F330w_AB', 'UVIS1_f336w_AB', 'WFC3_F336W_AB', 'u_SDSS_AB', 'U_Johnson', \
           'u3_CFHT_MC_AB', 'u1_CFHT_MC_AB', 'ACSWFC_F410w_AB', 'WFC3_FR388N_AB', 'ACSWFC_F435w_AB', 'WFC3_F438W_AB', 'UVIS1_f438w_AB', \
           'B3_Johnson', 'B2_Johnson', 'g_SDSS_AB', 'ACSWFC_F475w_AB', 'g3_CFHT_MC_AB', 'g1_CFHT_MC_AB', 'UVIS1_f555w_AB', \
           'WFC3_F555W_AB', 'ACSWFC_F555w_AB', 'UVIS1_f547m_AB', 'V_Johnson', 'ACSWFC_F606w_AB', 'UVIS1_f606w_AB', 'r_SDSS_AB', \
           'r1_CFHT_MC_AB', 'UVIS1_f625w_AB', 'ACSWFC_F625w_AB', 'r3_CFHT_MC_AB', 'UVIS1_f656n_AB', 'UVIS1_f657n_AB', 'UVIS1_f658n_AB', \
           'R_Cousins', 'R_Johnson', 'i_SDSS_AB', 'i2_CFHT_MC_AB', 'i3_CFHT_MC_AB', 'ACSWFC_F775w_AB', 'I_Cousins', 'UVIS1_f814w_AB', \
           'ACSWFC_F814w_AB', 'WFC3_F814W_AB', 'I_Johnson', 'z1_CFHT_MC_AB ', 'z_SDSS_AB', 'z3_CFHT_MC_AB', 'WFC3_F110W_AB', 'J_2Mass', \
           'J_Johnson', 'WFC3_F125W_AB', 'J_Palomar', 'WFC3_F160W_AB', 'H_Palomar', 'H_2MASS_AB', 'H_2Mass', 'Kprime_Cowie', 'Ks_CFHT_WC_AB', \
           'Ks_2Mass', 'K_Johnson', 'K_Palomar', 'L_Johnson', 'I3p6_IRAC', 'I4p5_IRAC', 'I5p7_IRAC', 'I7p9_IRAC', 'I12_IRAS', 'M24_MIPS', \
           'I25_IRAS', 'I60_IRAS', 'M70_MIPS', 'I100_IRAS', 'M160_MIPS', 'F(FUV)', 'F(NUV)', 'F(1500A)']
    bt.td3 = Table(names=s, dtype=['f4']*len(s))

    # Index file/table
    s   = ['log-age', 'CN_1', 'CN_2', 'Ca4227', 'G4300', 'Fe4383', 'Ca4455', 'Fe4531', 'Fe4668', 'Hbeta', 'Fe5015', \
           'Mg_1', 'Mg_2', 'Mg-b', 'Fe5270', 'Fe5335', 'Fe5406', 'Fe5709', 'Fe5782', 'Na-D', 'TiO_1', 'TiO_2', 'HdeltaA', 'HgammaA', \
           'HdeltaF ', 'HgammaF', 'B4_GC', 'B4_VN', 'CaII8498', 'CaII8542', 'CaII8662', 'MgI8807', 'H8_3889', 'H9_3835', 'H10_3798', \
           'BH-HK', 'BL1302', 'SiIV', 'BL1425', 'Fe1453', 'CIV1548a', 'CIV1548c', 'CIV1548e', 'BL1617', 'BL1664', 'BL1719', \
           'BL1853', 'FeII2402', 'BL2538', 'FeII2609', 'MgII', 'MgI', 'Mgwide', 'FeI', 'BL3096', 'CIVabs', 'HeIIems' ]
    bt.td4 = Table(names=s, dtype=['f4']*len(s))

    # Time step file/table
    if i == 0:
        # csp_galaxev output
        s   = ['age-yr', 'log_age-yr', 'MWA(yr)', 'uFWA(yr)', 'gFWA(yr)', 'rFWA(yr)', 'iFWA(yr)', 'zFWA(yr)', 'KFWA(yr)', 'MWLA(yr)', \
               'uFWLA(yr)', 'gFWLA(yr)', 'rFWLA(yr)', 'iFWLA(yr)', 'zFWLA(yr)', 'KFWLA(yr)', 'zu']
    else:
        # Other code output
        s   = ['age-yr', 'log_age-yr']
    bt.td5 = Table(names=s, dtype=['f8']*len(s))

    # Miscelaneous IMF information file/table
    s   = ['i', 'iseg', 'jvaz', 'jseg', 'x(i)', 'ml(i)', 'mu(i)', 'baux(i)', 'cn(i)', 'cc(i)']
    d   = ['i4', 'i4', 'i4', 'i4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4']
    bt.td6 = Table(names=s, dtype=d)

    # Information about bursts
    s   = ['Begin(yr)', 'Strength', 'File']
    d   = ['f4', 'f4', 'S']
    bt.td7 = Table(names=s, dtype=d)

    return kf,kp

def wrwa(t,wa,zu):
    # Writes flux weighted age to file open previously
    # First compute flux weighted age and log age
    for i in range(6):
        if wa[i+12] <= 0:
            wa[i]   = -99.
            wa[i+6] = -99.
        else:
            wa[i]   = wa[i]   / wa[i+12]
            wa[i+6] = wa[i+6] / wa[i+12]
    md = wa[20]
    if md <= 0:
        ma = -99.
        ml = -99.
    else:
        ma = wa[18]/md
        ml = wa[19]/md
    wa = [t,tlog(t),ma,wa[0],wa[1],wa[2],wa[3],wa[4],wa[5],ml,wa[6],wa[7],wa[8],wa[9],wa[10],wa[11],zu]
    return wa

def tlog(t):
    # avoids log of 0
    if t <= 0:
        tlog = 0.
    else:
        tlog = np.log10(t)
    return tlog

def rf_props(t,x,y,b,sf,k):
    # Gets rest frame properties of galaxy model at log(age) = tl
    # Properties are computed and saved in the order used in the model fits file
    global bolflux,wa

    # Store time step and logarithmic time step
    tl = tlog(t)
    if k==0:
        # csp_galaxev, store flux weighted age data as well
        wt = wrwa(t,bt.wa,bt.zu)
        bt.td5.add_row(wt)
    else:
        # other codes
        bt.td5.add_row((t,tl))

    # Compute absolute magnitude in filters in array kf at z = 0
    bt.bolflux = b[0]
   #bolflux = np.trapz(y,x)
    vm = fl.absmag(x,y,bt.kf)
    vm.insert(0,tl)
    bt.td3.add_row(vm)

    # Compute Lick+ line strength indices
    xi = ix.indx(x,y,57)
    xi.insert(0,tl)
    bt.td4.add_row(xi)

    def insert(p,a,v):
        # Inserts at position p of array a values in array v
        for i in range(len(v)):
            a.insert(p,v[i])
        return a

    # Compute physical properties and store them in array pq
    #   Various rates, number of remnants, mass of different galaxy components
    pq=b
    #   Evolutionary flux
    evflux = pq[len(pq)-5]              # Evolutionary flux
    if bt.bolflux > 0:
        evf = evflux/bt.bolflux      # Evolutionary flux/ Bolometric flux
    else:
        evf = 0.
    pq.insert(16,evf)
    #   Mass-to-light ratios
    mstr = pq[10] ; mtot = mstr + pq[11]
    bmag = vm[21] ; vmag = vm[31] ; kmag = vm[69]
    mlb1,mlv1,mlk1  = fl.mass2light(mstr,bmag,vmag,kmag)
    mlb, mlv, mlk   = fl.mass2light(mtot,bmag,vmag,kmag)
    pq = insert(15,pq,[mlk1,mlv1,mlb1,mlk,mlv,mlb,mtot])
    #   X-ray luminosity
    pq.insert(10,fl.lx(x,y))
    #   Ionizing photons
    io = fl.ionphot(x,y) ; io.insert(0,tl)
    #   Spectral breaks
    br = ix.breaks(x,y,['912','vn','sdss','b83'])
    #   Concatenate arrays
    del pq[0]
    r2 = np.concatenate((io,br,pq))
    #    SFR
    r2[23] = sf
    bt.td2.add_row(r2)

def percent(i,n,label,p):
    # Reports percentage done
    if p<0:
        print()
        print(label+':')
        print(' % done: ...10...20...30...40...50...60...70...80...90..100 (' + str(n) + ' steps)')
        p=0
    qr=float(i)/float(n)*50.
    qr=int(qr)
    if qr>p:
        l=''
        for j in range(qr+1):
            l=l+'.'
        print('         '+l, end='\r')
        p=qr
    if i>=n:
        print()
    return p


import numpy as np
import basics.basics as bs

def indx(x,y,n):
    # Compute line strength indices

    def bh_index(j,x,y):
        # Computes value of Brodie and Hanes (1986) index J directly from sed
        # Follows procedure outlined by Brodie and Hanes (1986, ApJ 300, 258)
        # See their eqs. 1-2, and Table 2 in pag. 260.

        # Define wavelength range for selected index
        if j==1:
            w = [ 3800., 4000., 4000., 4200. ] # idxlbl = 'Delta (UV Blanketing)'
        elif j==2:
            w = [ 3785., 3795., 3810., 3825. ] # idxlbl = 'H10 (Balmer line)'
        elif j==3:
            w = [ 3810., 3825., 3850., 3865. ] # idxlbl = 'H9 (Balmer line)'
        elif j==4:
            w = [ 3865., 3885., 3910., 3925. ] # idxlbl = 'Hpsi (Balmer line)'
        elif j==5:
            w = [ 3785., 3810., 3910., 3925. ] # idxlbl = 'CN (UV cyanogen)'
        elif j==6:
            w = [ 3910., 3925., 3995., 4015. ] # idxlbl = 'HK index (Ca II H and K lines)'
        elif j==7:
            w = [ 4075., 4090., 4125., 4140. ] # idxlbl = 'Hdelta (Balmer line)'
        elif j==8:
            w = [ 4200., 4215., 4245., 4260. ] # idxlbl = 'Ca (Ca I line)'
        elif j==9:
            w = [ 4275., 4285., 4315., 4325. ] # idxlbl = 'G (CH G band)'
        elif j==10:
            w = [ 4315., 4325., 4360., 4375. ] # idxlbl = 'Hgamma (Balmer line)'
        elif j==11:
            w = [ 4800., 4830., 4890., 4920. ] # idxlbl = 'Hbeta (Balmer line)'
        elif j==12:
            w = [ 5125., 5150., 5195., 5220. ] # idxlbl = 'Mg (Mg b triplet)'
        elif j==13:
            w = [ 4740., 4940., 5350., 5550. ] # idxlbl = 'MH (Mg b + MgH)'
        elif j==14:
            w = [ 5225., 5250., 5280., 5305. ] # idxlbl = 'FC (Fe I + Ca I)'
        elif j==15:
            w = [ 5385., 5865., 5920., 5950. ] # idxlbl = 'Na (Na D lines)'
        elif j==16:
            w = [ 5950., 6100., 6330., 6480. ] # idxlbl = 'FTC (Fe I + TiO + Ca I)'
        elif j==17:
            w = [ 6480., 6615., 6575., 6610. ] # idxlbl = 'Halpha (Balmer line)'
    
        # Compute index according to BH definition
        # Mean flux through blue window
        fb=bs.trapz2(x,y,w[0],w[1])/(w[1]-w[0])		# fb=trapz3(x,f,n,w(1),w(2),ierr)/(w(2)-w(1))
        # Mean flux through red window
        fr=bs.trapz2(x,y,w[2],w[3])/(w[3]-w[2])		# fr=trapz3(x,f,n,w(3),w(4),ierr)/(w(4)-w(3))
        # Compute index
        if j==1:
            if fb <= 0:
                bh_index = -99.
            else:
                bh_index = 2.5*np.log10(fr/fb)
        else:
            # Mean flux through central window
            fc=bs.trapz2(x,y,w[1],w[2])/(w[2]-w[1])		# fc=trapz3(x,f,n,w(2),w(3),ierr)/(w(3)-w(2))
            if fb+fr <= 0:
                bh_index = -99
            else:
                bh_index = 1.0 - 2.0*fc/(fb+fr)
        return bh_index

    def lsindx_sed(j,x,y):

        def fc(v):
            # Equation of straight line joining (wb,fb) and (wr,fr)
            fc = ( (wr-v)*fb + (v-wb)*fr ) / (wr-wb)
            return fc

        if j==1:
            im = 1 ; w = [ 4080.125, 4117.625, 4244.125, 4284.125, 4142.125, 4177.125 ] ; idxlbl = 'Lick index 01: CN 1 molecular band (mag)'
        elif j==2:
            im = 1 ; w = [ 4083.875, 4096.375, 4244.125, 4284.125, 4142.125, 4177.125 ] ; idxlbl = 'Lick index 02: CN 2 molecular band (mag)'
        elif j==3:
            im = 0 ; w = [ 4211.000, 4219.750, 4241.000, 4251.000, 4222.250, 4234.750 ] ; idxlbl = 'Lick index 03: Ca(4227 index (Ang)'
        elif j==4:
            im = 0 ; w = [ 4266.375, 4282.625, 4318.875, 4335.125, 4281.375, 4316.375 ] ; idxlbl = 'Lick index 04: G4300 index (Ang)'
        elif j==5:
            im = 0 ; w = [ 4359.125, 4370.375, 4442.875, 4455.375, 4369.125, 4420.375 ] ; idxlbl = 'Lick index 05: Fe4383 index (Ang)'
        elif j==6:
            im = 0 ; w = [ 4445.875, 4454.625, 4477.125, 4492.125, 4452.125, 4474.625 ] ; idxlbl = 'Lick index 06: Ca4455 index (Ang)'
        elif j==7:
            im = 0 ; w = [ 4504.250, 4514.250, 4560.500, 4579.250, 4514.250, 4559.250 ] ; idxlbl = 'Lick index 07: Fe4531 index (Ang)'
        elif j==8:
            im = 0 ; w = [ 4611.500, 4630.250, 4742.750, 4756.500, 4634.000, 4720.250 ] ; idxlbl = 'Lick index 08: Fe4668 index (Ang)'
        elif j==9:
            im = 0 ; w = [ 4827.875, 4847.875, 4876.625, 4891.625, 4847.875, 4876.625 ] ; idxlbl = 'Lick index 09: H beta index (Ang)'
        elif j==10:
            im = 0 ; w = [ 4946.500, 4977.750, 5054.000, 5065.250, 4977.750, 5054.000 ] ; idxlbl = 'Lick index 10: Fe5015 index (Ang)'
        elif j==11:
            im = 1 ; w = [ 4895.125, 4957.625, 5301.125, 5366.125, 5069.125, 5134.125 ] ; idxlbl = 'Lick index 11: Mg 1 index (mag)'
        elif j==12:
            im = 1 ; w = [ 4895.125, 4957.625, 5301.125, 5366.125, 5154.125, 5196.625 ] ; idxlbl = 'Lick index 12: Mg 2 index (mag)'
        elif j==13:
            im = 0 ; w = [ 5142.625, 5161.375, 5191.375, 5206.375, 5160.125, 5192.625 ] ; idxlbl = 'Lick index 13: Mg b index (Ang)'
        elif j==14:
            im = 0 ; w = [ 5233.150, 5248.150, 5285.650, 5318.150, 5245.650, 5285.650 ] ; idxlbl = 'Lick index 14: Fe5270 index (Ang)'
        elif j==15:
            im = 0 ; w = [ 5304.625, 5315.875, 5353.375, 5363.375, 5312.125, 5352.125 ] ; idxlbl = 'Index 15: Fe5335 index (Ang)'
        elif j==16:
            im = 0 ; w = [ 5376.250, 5387.500, 5415.000, 5425.000, 5387.500, 5415.000 ] ; idxlbl = 'Lick index 16: Fe5406 index (Ang)'
        elif j==17:
            im = 0 ; w = [ 5672.875, 5696.625, 5722.875, 5736.625, 5696.625, 5720.375 ] ; idxlbl = 'Lick index 17: Fe5709 index (Ang)'
        elif j==18:
            im = 0 ; w = [ 5765.375, 5775.375, 5797.875, 5811.625, 5776.625, 5796.625 ] ; idxlbl = 'Lick index 18: Fe5782 index (Ang)'
        elif j==19:
            im = 0 ; w = [ 5860.625, 5875.625, 5922.125, 5948.125, 5876.875, 5909.375 ] ; idxlbl = 'Lick index 19: Na D index (Ang)'
        elif j==20:
            im = 1 ; w = [ 5816.625, 5849.125, 6038.625, 6103.625, 5936.625, 5994.125 ] ; idxlbl = 'Lick index 20: TiO 1 index (mag)'
        elif j==21:
            im = 1 ; w = [ 6066.625, 6141.625, 6372.625, 6415.125, 6189.625, 6272.125 ] ; idxlbl = 'Lick index 21: TiO 2 index (mag)'
        elif j==22:
            im = 0 ; w = [ 4041.600, 4079.750, 4128.500, 4161.000, 4083.500, 4122.250 ] ; idxlbl = 'Hdelta_A (Ang)' # Index defined by Worthey and Ottaviani (1997)
        elif j==23:
            im = 0 ; w = [ 4283.500, 4319.750, 4367.250, 4419.750, 4319.750, 4363.500 ] ; idxlbl = 'Hgamma_A (Ang)' # Index defined by Worthey and Ottaviani (1997)
        elif j==24:
            im = 0 ; w = [ 4057.250, 4088.500, 4114.750, 4137.250, 4091.000, 4112.250 ] ; idxlbl = 'Hdelta_F (Ang)' # Index defined by Worthey and Ottaviani (1997)
        elif j==25:
            im = 0 ; w = [ 4283.500, 4319.750, 4354.750, 4384.750, 4331.250, 4352.250 ] ; idxlbl = 'Hgamma_F (Ang)' # Index defined by Worthey and Ottaviani (1997)
        elif j==26:		# No. 30 in fortran version
            return sbreak(x,y,'gc')											# Gorgas and Cardiel D4000 index
        elif j==27:		# No. 31 in fortran version
            return sbreak(x,y,'vn')											# B4000VN index
        elif j==28:		# No. 26 in fortran version
            im = 0 ; w = [ 8447.500, 8462.500, 8842.500, 8857.500, 8483.000, 8513.000 ] ; idxlbl = 'CaII8498 (Ang)' # Index Ca1 defined by Diaz, Terlevich, and Terlevich (DTT, 1989)
        elif j==29:		# No. 27 in fortran version
            im = 0 ; w = [ 8447.500, 8462.500, 8842.500, 8857.500, 8527.000, 8557.000 ] ; idxlbl = 'CaII8542 (Ang)' # Index Ca2 defined by Diaz, Terlevich, and Terlevich (DTT, 1989)
        elif j==30:		# No. 28 in fortran version
            im = 0 ; w = [ 8447.500, 8462.500, 8842.500, 8857.500, 8647.000, 8677.000 ] ; idxlbl = 'CaII8662 (Ang)' # Index Ca3 defined by Diaz, Terlevich, and Terlevich (DTT, 1989)
        elif j==31:		# No. 29 in fortran version
            im = 0 ; w = [ 8775.000, 8787.000, 8845.000, 8855.000, 8799.500, 8814.500 ] ; idxlbl = 'MgI8807  (Ang)' # Index MgI defined by Diaz, Terlevich, and Terlevich (DTT, 1989)
        elif j==32:
            im = 0 ; w = [ 3855.000, 3865.000, 3905.000, 3915.000, 3870.000, 3900.000 ] ; idxlbl = 'H8_3889  (Ang)' # Index defined by Delphine Marcillac (2004)
        elif j==33:
            im = 0 ; w = [ 3810.000, 3820.000, 3855.000, 3865.000, 3825.000, 3845.000 ] ; idxlbl = 'H9_3835  (Ang)' # Index defined by Delphine Marcillac (2004)
        elif j==34:
            im = 0 ; w = [ 3775.000, 3785.000, 3810.000, 3820.000, 3790.000, 3805.000 ] ; idxlbl = 'H10_3798 (Ang)' # Index defined by Delphine Marcillac (2004)
        elif j==35:
            return bh_index(6,x,y)											# Brodie and Hanes HK index
        elif j==36:		# Fanelli et al. UV indices 1-21 start here:
            im = 0 ; w = [ 1270.000, 1290.000, 1345.000, 1365.000, 1292.000, 1312.000 ] ; idxlbl = 'Index 01: BL 1302 (A)'                     # Si III, Si II, O I
        elif j==37:
            im = 0 ; w = [ 1345.000, 1365.000, 1475.000, 1495.000, 1387.000, 1407.000 ] ; idxlbl = 'Index 02: Si IV (A)'                       # Si IV 1393.8; 1402.8
        elif j==38:
            im = 0 ; w = [ 1345.000, 1365.000, 1475.000, 1495.000, 1415.000, 1435.000 ] ; idxlbl = 'Index 03: BL 1425 (A)'                     # C II 1429, Si III 1417, Fe IV, Fe V
        elif j==39:
            im = 0 ; w = [ 1345.000, 1365.000, 1475.000, 1495.000, 1440.000, 1466.000 ] ; idxlbl = 'Index 04: Fe 1453 (A)'                     # Fe V + 20 additional lines
        elif j==40:
            im = 0 ; w = [ 1500.000, 1520.000, 1577.000, 1597.000, 1530.000, 1550.000 ] ; idxlbl = 'Index 05: C IV 1548 absorption (A)'        # C IV 1548, in absorption
        elif j==41:
            im = 0 ; w = [ 1500.000, 1520.000, 1577.000, 1597.000, 1540.000, 1560.000 ] ; idxlbl = 'Index 06: C IV 1548 central band (A)'      # C IV 1548, central band
        elif j==42:
            im = 0 ; w = [ 1500.000, 1520.000, 1577.000, 1597.000, 1550.000, 1570.000 ] ; idxlbl = 'Index 07: C IV 1548 emission (A)'          # C IV 1548, in emission
        elif j==43:
            im = 0 ; w = [ 1577.000, 1597.000, 1685.000, 1705.000, 1604.000, 1630.000 ] ; idxlbl = 'Index 08: BL 1617 (A)'                     # Fe IV
        elif j==44:
            im = 0 ; w = [ 1577.000, 1597.000, 1685.000, 1705.000, 1651.000, 1677.000 ] ; idxlbl = 'Index 09: BL 1664 (A)'                     # C I 1656.9, Al II 1670.8
        elif j==45:
            im = 0 ; w = [ 1685.000, 1705.000, 1803.000, 1823.000, 1709.000, 1729.000 ] ; idxlbl = 'Index 10: BL 1719 (A)'                     # N IV 1718.6, Si IV 1722.5,1727.4, Al II
        elif j==46:
            im = 0 ; w = [ 1803.000, 1823.000, 1885.000, 1915.000, 1838.000, 1868.000 ] ; idxlbl = 'Index 11: BL 1853 (A)'                     # Al II, Al III, Fe II, Fe III
        elif j==47:
            im = 0 ; w = [ 2285.000, 2325.000, 2432.000, 2458.000, 2382.000, 2422.000 ] ; idxlbl = 'Index 12: Fe II 2402 (A)'                  # Fe II 2402
        elif j==48:
            im = 0 ; w = [ 2432.000, 2458.000, 2562.000, 2588.000, 2520.000, 2556.000 ] ; idxlbl = 'Index 13: BL 2538 (A)'                     # Perhaps Fe I
        elif j==49:
            im = 0 ; w = [ 2562.000, 2588.000, 2647.000, 2673.000, 2596.000, 2622.000 ] ; idxlbl = 'Index 14: Fe II 2609 (A)'                  # Fe II 2609
        elif j==50:
            im = 0 ; w = [ 2762.000, 2782.000, 2818.000, 2838.000, 2784.000, 2814.000 ] ; idxlbl = 'Index 15: Mg II (A)'
        elif j==51:
            im = 0 ; w = [ 2818.000, 2838.000, 2906.000, 2936.000, 2839.000, 2865.000 ] ; idxlbl = 'Index 16: Mg I (A)'
        elif j==52:
            im = 0 ; w = [ 2470.000, 2670.000, 2930.000, 3130.000, 2670.000, 2870.000 ] ; idxlbl = 'Index 17: Mg wide (A)'
        elif j==53:
            im = 0 ; w = [ 2906.000, 2936.000, 3031.000, 3051.000, 2965.000, 3025.000 ] ; idxlbl = 'Index 18: Fe I (A)'
        elif j==54:
            im = 0 ; w = [ 3031.000, 3051.000, 3115.000, 3155.000, 3086.000, 3106.000 ] ; idxlbl = 'Index 19: BL 3096 (A)'                     # Al I 3092, Fe I 3091.6
        elif j==55:
            im = 0 ; w = [ 1510.000, 1525.000, 1560.000, 1575.000, 1525.000, 1548.000 ] ; idxlbl = 'CIV stellar absorption (A)'
        elif j==56:
            im = 0 ; w = [ 1600.000, 1650.000, 1650.000, 1700.000, 1635.000, 1645.000 ] ; idxlbl = 'HeII stellar emission (A)'
    
        # Compute mean height inside blue pseudo continuum band
        fb=bs.trapz2(x,y,w[0],w[1])/(w[1]-w[0])
        wb=(w[0]+w[1])/2.

        # Compute mean height inside red pseudo continuum band
        fr=bs.trapz2(x,y,w[2],w[3])/(w[3]-w[2])
        wr=(w[2]+w[3])/2.

        # Check if zero flux
        if fb == 0 and fr == 0:
            gw_ix_sed = -99.
            return gw_ix_sed

        # Compute equivalent width according to eq. (2) of Trager at el. (1998, ApJS 116, 1)
        # Locate central band in array x: w(5) and w(6)
        i1, i2 = np.searchsorted(x,[w[4],w[5]])
        # Fill in auxiliary arrays from i1-5 to i2+5 with ratio f/fc
        u=[]
        v=[]
        for i in range (i1-5,i2+6):
            u.append(x[i])
            v.append(y[i]/fc(x[i]))
        # Compute integral of ratio f(i)/fc(x(i)) from w(5) to w(6)
        fm=bs.trapz2(u,v,w[4],w[5])
        # Compute index
        if im == 0:
            gw_ix_sed = w[5] - w[4] - fm
        else:
            fm=fm / (w[5] - w[4])
            gw_ix_sed = -2.5*np.log10(fm)

        # Return indez value
        return gw_ix_sed

    ix=[]
    for k in range(1,n):
        ix.append(lsindx_sed(k,x,y))
    return ix

def sbreak(x,y,c):

    # Measures various breaks in sed in (x,y) according to choice in 'c'
    # The break is defined as the ratio of the average Fnu flux densities in the ranges:
    # (w[2]:w[3])/(w[0]:w[1])

    # Select wavelength range
    if c=='b83':
        # 4000 A break
        # (4050-4250)/(3750-3950) (Bruzual (1983), Hamilton 1985 ApJ 297, 371).
        w = [ 3750., 3950., 4050., 4250. ]
    elif c=='vn':
        # 4000 A break
        # Very narrow definition, via S Charlot from MPA
        w = [ 3850., 3950., 4000., 4100. ]
    elif c=='sdss':
        # 4000 A break
    	# Uses definition by Stoughton et al. (Ap. J), for SDSS data
        w = [ 3751., 3951., 4051., 4251. ]
    elif c=='gc':
        # 4000 A break
        # J. Gorgas and N. Cardiel definition
        w = [ 3750., 3950., 4050., 4250. ]
    elif c=='912':
        # 912 A break
        # (800-900)/(1000-1100) (Bruzual (1983).
        w = [ 800., 900., 1000., 1100. ]
    i1, i2, i3, i4 = np.searchsorted(x,[w[0],w[1],w[2],w[3]])

    # Compute average flux below in Fnu units
    v = []
    f = []
    for i in range(i1,i2+1):
        v.append(x[i])
        f.append(y[i]*x[i]**2)
    fl = np.trapz(f,v)/(w[1]-w[0])

    # Compute average flux above in Fnu units
    v = []
    f = []
    for i in range(i3,i4+1):
        v.append(x[i])
        f.append(y[i]*x[i]**2)
    fr = np.trapz(f,v)/(w[3]-w[2])

    # Compute break amplitude
    if fl <= 0:
        sbreak = -99.
    else:
        sbreak=fr/fl
    return sbreak

def breaks(x,y,b):
    # computes spectral breaks listed in b
    br=[]
    for i in range(len(b)):
        br.append(sbreak(x,y,b[i]))
    return br


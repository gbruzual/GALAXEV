import os
import sys
import numpy         as np
import builtins      as bt
import basics.basics as bs
import bcfits.bcfits as bc
import bchead.bchead as hd
import common.common as cn
import bcfilt.bcfilt as fl
import cosmol.cosmol as cl

def updp(i,file1):
    # define cosmological parameters
    global h, omega, omega_lambda, clambda, q, tu, ttg, zf

    if i == 0:
        # Recover default values of cosmological parameters
        h, q, clambda, tu, ttg, zf, omega, omega_lambda = bt.cosmol
        # Return cosmological parameters
        # return h, q, clambda, tu, ttg, zf, omega, omega_lambda

    elif i == 1:
        print()
        s = input(' BC_GALAXEV model sed in file [' + bs.lnam(file1) + '] = ')
        if len(s) <= 0:
            s = bs.lnam(file1)
        s = bs.fcheck(s)
        return s

    elif i == 2:
        # Change default value of cosmological parameters
        print()
        print(' COSMOLOGY: Modify values of cosmological parameters if desired (separated by commas):')
        s = input('    Enter Ho [' + str('{:.1f}'.format(h)) + '], Omega [' + str('{:.3f}'.format(omega)) + '], Omega_lambda [' + str('{:.3f}'.format(omega_lambda)) +'] = ')
        if (len(s) > 0):
            s = s.replace(" ", ",")
            s = s.split(",")
            for i in range(4-len(s)):
                s.append("0.")
            s[0], s[1], s[2] = float(s[0]), float(s[1]), float(s[2])
            if s[0] > 0.:
                h = s[0]
            if s[1] > 0.:
                omega = s[1]
            if s[2] > 0.:
                omega_lambda = s[2]
            s = str('{:.2f}'.format(h)) + ',' + str('{:.4f}'.format(omega)) + ',' + str('{:.4f}'.format(omega_lambda))
            clambda,q = cl.cosmol_c(h,omega,omega_lambda)
            tu=cl.tuniverse(h,q,0.,clambda)
            ttg=min(tu,ttg)
            zf=cl.zx(ttg,h,q,clambda)
        return s

    elif i == 3:
        # Change age of galaxies
        clambda,q = cl.cosmol_c(h,omega,omega_lambda)
        tu=cl.tuniverse(h,q,0.,clambda)
       #print()
        print('    Age of this universe = tu = ' + str('{:.2f}'.format(tu)) +' Gyr')
        s = input('      Enter age of galaxy today = tg [' + str('{:.1f}'.format(ttg)) + ' Gyr] = ')
        if (len(s) > 0):
            t = float(s)
            if (t > 0.):
                ttg = t
        ttg=min(tu,ttg)
        zf=cl.zx(ttg,h,q,clambda)
        return s

    elif i == 4:
       #print()
        print('    Cosmological model in use:')
        print('      Ho = ' + str('{:.1f}'.format(h)) + '   Omega = ' + str('{:.2f}'.format(omega)) + '   Omega_lambda = ' + str('{:.2f}'.format(omega_lambda)) +
                 '   qo = ' + str('{:.4f}'.format(q)) + '   Lambda = ' + str('{:.3E}'.format(clambda)) + '  tu = ' + str('{:.2f}'.format(tu)) + ' Gyr   tg = ' + 
                     str('{:.2f}'.format(ttg)) + ' Gyr    zf = ' + str('{:.2f}'.format(zf)))
       #print('      tu = ' + str('{:.2f}'.format(tu)) + ' Gyr   tg = ' + str('{:.2f}'.format(ttg)) + ' Gyr    zf = ' + str('{:.2f}'.format(zf)))
        return

def header3(f,e):
    # Writes header in BC tables
    a = '#          '
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
            if 'LISTED' in s:
                s = 'I      LISTED: Observer frame properties for the indicated cosmology                                                     I'
            f.write(a + s + '\n')
            i=i+2
    f.write('#' + '\n')

def zage(t,ttg):
    global kl, zz
    h, q, clambda, *nouse = bt.cosmol
    # Get z corresponding to each galaxy age
   #kl = np.searchsorted(t, ttg*1.E9) - 1
    kl = np.searchsorted(t, ttg*1.E9)
    tg = t[kl]
    zz=[]
    for k in range(kl+1):
        i1=kl-k
        t1=(tg-t[i1])*1.e-9
        z1=cl.zx(t1,h,q,clambda)
        zz.append(z1)
    zz=np.flip(zz)
    return zz

def zsed(z,t,f):
    # Returns SED corresponding to redshift z according to array zz
    # zz contains the redshift corresponding to a given look-back-time tx
    h, q, clambda, *nouse = bt.cosmol
    if z >= zz[0]:
        # Highest redshift, use first time step sed
        # ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # ! If z > zz(i1), use first sed in model (=> no evolution at this z)
        # ! Modified as suggested by Rogier Windhorst to allow high z
        # ! colors to be computed. This means that after the t=0
        # ! sed from galaxev is used, only the k-correction will be
        # ! computed (no more evolution possible because galaxy
        # ! did not exist at higher z''s)
        # ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        yz=bc.ft(t[0]*1.E-9,f,t)
        tx = ttg
    else:
        tg = t[kl]
        for k in range(kl):
            i1=kl-k
            if zz[i1] <= z and z < zz[i1-1]:
                i2=i1-1
                t1=(tg-t[i1])*1.e-9
                t2=(tg-t[i2])*1.e-9
                # Look-back-time for redshift z
                tu=cl.tuniverse(h,q,0.,clambda)
                tx=tu - cl.tuniverse(h,q,z,clambda)
                a1=(t2-tx)/(t2-t1)
                a2=1.-a1
                y1=bc.ft(t[i1],f,t)
                y2=bc.ft(t[i2],f,t)
                yz=a1*y1+a2*y2
                break
    return tx, yz

def zk(kz,zf):
    #! Returns redshift used at step kz
    if kz <= 21:
        zk =  0.000 + (kz  -1)*0.005
    elif kz <= 97:
        zk =  0.100 + (kz -21)*0.025
    elif kz <= 177:
        zk =  2.000 + (kz -97)*0.100
    else:
       #zk = 10.000 + (kz-177)*0.020
        zk = 10.000 + (kz-177)*0.200
    if zk > zf:
        zk=zf
    return zk

def cmev(file1,jupy):
    # python version of cm_evolution.f code (alias cmev)
    global h, omega, omega_lambda, clambda, q, tu, ttg, zf

    def meaning(o,nf):
        global n1,n2

        # Prints meaning of different quantities in output files
        if o==o1 or o==o2:
            l='# Cosmology: Ho = '+f'{h:4.0f}'+', Omega ='+f'{omega:5.2f}'+', Omega_lambda ='+f'{omega_lambda:5.2f}'+', qo = '+f'{q:7.4f}'+', Lambda ='+f'{clambda:10.3E}'+', tu = '+f'{tu:5.2f} Gyr'+', tg = '+f'{ttg:5.2f} Gyr'+',  zf = '+f'{zf:6.2f}\n'
            o.write(l)
            o.write('#\n')
            o.write('# Meaning of printed variables (see Tables 5 and 6 of bc03.pdf for details):\n')
            o.write('#   z        = redshift (z=0 entry in table corresponds to a distance of 10 pc)\n')
            o.write('#   ltt      = light travel time from z=0 to z in Gyr\n')
            o.write('#   tz       = age of galaxy at redshit = z\n')
            o.write('#   dm       = cosmological distance modulus in magnitude units\n')
            o.write('#   M_rf     = M_rf(z) in filter ' + nf + ' (Vega system)\n')
            o.write('#   M_ne     = M_ne(z) in filter ' + nf + ' (Vega system)\n')
            o.write('#   M_ev     = M_ev(z) in filter ' + nf + ' (Vega system)\n')
            o.write('#   m_ev     = apparent magnitude M_ev(z) + dm(z) in filter ' + nf + ' (Vega system)\n')
            o.write('#   M_AB_rf  = AB M_rf(z) in filter ' + nf + ' (AB system)\n')
            o.write('#   M_AB_ne  = AB M_ne(z) in filter ' + nf + ' (AB system)\n')
            o.write('#   M_AB_ev  = AB M_ev(z) in filter ' + nf + ' (AB system)\n')
            o.write('#   m_AB_ev  = apparent AB magnitude M_ev(z) + dm(z) in filter ' + nf + ' (AB system)\n')
            o.write('#   e+k_cor  = (e+k)-correction in filter '+ nf + '\n')
            o.write('#   k_cor_ev = k-correction computed for sed at redshift z in filter '+ nf + '\n')
            o.write('#   k_cor_ne = k-correction computed for sed at redshift 0 in filter '+ nf + '\n')
            o.write('#\n')
            if o==o1:
                n1=int(nf)
                l='# Filter:                      '+f'{n1:9d}{n1:9d}{n1:9d}{n1:9d}   {n1:9d}{n1:9d}{n1:9d}{n1:9d}   {n1:9d}{n1:9d}{n1:9d}\n'
                o.write(l)
            else:
                n2=int(nf)
                l='# Filter:                      '+f'{n2:9d}{n2:9d}{n2:9d}{n2:9d}   {n2:9d}{n2:9d}{n2:9d}{n2:9d}   {n2:9d}{n2:9d}{n2:9d}\n'
                o.write(l)
            o.write('#  z      LTT     tz      dm         M_rf     M_ne     M_ev     m_ev      M_AB_rf  M_AB_ne  M_AB_ev  m_AB_ev     e+k_cor  k_cor_ev k_cor_ne\n')
            o.write('#         Gyr     Gyr     mag        mag      mag      mag      mag        AB_mag   AB_mag   AB_mag   AB_mag       mag      mag      mag\n')
            o.write('# (1)     (2)     (3)     (4)        (5)      (6)      (7)      (8)         (9)      (10)     (11)     (12)        (13)     (14)     (15)\n')
        elif o==o3:
            l='# Cosmology: Ho = '+f'{h:4.0f}'+', Omega ='+f'{omega:5.2f}'+', Omega_lambda ='+f'{omega_lambda:5.2f}'+', qo = '+f'{q:7.4f}'+', Lambda ='+f'{clambda:10.3E}'+', tu = '+f'{tu:5.2f} Gyr'+', tg = '+f'{ttg:5.2f} Gyr'+',  zf = '+f'{zf:6.2f}\n'
            o.write(l)
            o.write('#\n')
            o.write('# Meaning of printed variables (see Tables 5 and 6 of bc03.pdf for details):\n')
            o.write('#   z        = redshift (z=0 entry in table corresponds to a distance of 10 pc)\n')
            o.write('#   ltt      = light travel time from z=0 to z in Gyr\n')
            o.write('#   tz       = age of galaxy at redshit = z\n')
            o.write('#   dm       = cosmological distance modulus in magnitude units\n')
            o.write('#   C_rf     = color M1_rf(z) - M2_rf(z) (Vega system)\n')
            o.write('#   C_ne     = color M1_ne(z) - M2_ne(z) (Vega system)\n')
            o.write('#   C_ev     = color M1_ev(z) - M2_ev(z) (Vega system)\n')
            o.write('#   C_AB_rf  = AB color M1_rf(z) - M2_rf(z) (AB system)\n')
            o.write('#   C_AB_ne  = AB color M1_ne(z) - M2_ne(z) (AB system)\n')
            o.write('#   C_AB_ev  = AB color M1_ev(z) - M2_ev(z) (AB system)\n')
            o.write('#   e+k_cor  = color (e+k)-correction\n')
            o.write('#   k_cor_ev = color k-correction computed for sed at redshift z\n')
            o.write('#   k_cor_ne = color k-correction computed for sed at redshift 0\n')
            o.write('#\n')
            l='# Filters:                       '+f'{n1:5d}-{n2:3d}{n1:5d}-{n2:3d}{n1:5d}-{n2:3d}   {n1:5d}-{n2:3d}{n1:5d}-{n2:3d}{n1:5d}-{n2:3d}   {n1:5d}-{n2:3d}{n1:5d}-{n2:3d}{n1:5d}-{n2:3d}\n'
            o.write(l)
            o.write('#  z      LTT     tz      dm         C_rf     C_ne     C_ev      C_AB_rf  C_AB_ne  C_AB_ev     e+k_cor  k_cor_ev k_cor_ne\n')
            o.write('#         Gyr     Gyr     mag        mag      mag      mag        AB_mag   AB_mag   AB_mag       mag      mag      mag\n')
            o.write('# (1)     (2)     (3)     (4)        (5)      (6)      (7)         (8)      (9)      (10)        (11)     (12)     (13)\n')

    global f, t, fd, iread

    # Init variables
    bt.iread = True     # Read filter file only on first call
    p=-1                # Needed by 'percent' function on first call

    if not jupy:
        # Ask for input parameters to run cm_evolution program
        bs.bchead()
        f1 = input(' BC_GALAXEV SSP sed in file [' + bs.lnam(file1) + '] = ')
        if len(f1) <= 0:
            file1 = bs.lnam(file1)
        else:
            file1 = bs.zrep(file1,f1)

        # Change default cosmology
        updp(0,file1)
        s = updp(2,file1)
        print(s)
        s = updp(3,file1)
        print(s)
        updp(4,file1)
        #s = updp(1,file1)
        #print(s)
        print()
        print()
        print(' PHOTOMETRY: Enter up to 2 filters separated by commas.')
        af = input('    Compute magnitude/color through filter numbers [14,15] = ')
        if af=='':
            af = '14,15'
    else:
        # Recover default values of cosmological parameters
        h, q, clambda, tu, ttg, zf, omega, omega_lambda = bt.cosmol

        # Get parameters from data entered in widget
        # ['cb2019_z017_chab_hr_xmilesi_ssp.fits', '71,0.27,0.73', '13.5', '20,21']
        a = file1
        a[1] = a[1].split(',')
        h = float(a[1][0])
        omega = float(a[1][1])
        omega_lambda = float(a[1][2])
        clambda,q = cl.cosmol_c(h,omega,omega_lambda)
        tu=cl.tuniverse(h,q,0.,clambda)
        ttg=float(a[2])
        ttg=min(tu,ttg)
        zf=cl.zx(ttg,h,q,clambda)
        af=a[3]
        file1=a[0]

    # Read BC/CB model in fits file
    w,f,t,e,*nouse = bc.bcfits(file1)
    hd.multhead(0,file1,t,w,e)

    # Filter pair
    af = af.split(',')
    af = np.array(af)
    kf = [int(x) for x in af]

    # Build output file name and write header into it
    for i in range(len(kf)):
        ext='F000'
        mf=af[i]
        jf=int(mf)
        if jf < 10:
            ext = 'F00' + mf
        elif jf < 100:
            ext = 'F0' + mf
        else:
            ext = 'F' + mf
        oo = os.environ.get('glxout') + '/'
        oo = oo.replace('//','/')
        if i==0:
            imag=1
            ext1=ext
            ext='magnitude_' + ext1
            name=file1.replace('fits',ext)
            o1 = open(oo + name,'w')
            header3(o1,e)
            meaning(o1,mf)
            h1 = ' Output written to file(s): ' + oo + name
        elif mf !='0':
            imag=2
            ext2=ext
            ext='magnitude_' + ext2
            name=file1.replace('fits',ext)
            o2 = open(oo + name,'w')
            header3(o2,e)
            meaning(o2,mf)
            h2 = '                            ' + oo + name
            ext='color_' + ext1 + '_' + ext2
            name=file1.replace('fits',ext)
            o3 = open(oo + name,'w')
            header3(o3,e)
            meaning(o3,mf)
            h3 = '                            ' + oo + name

    # Get z corresponding to each look-back-time in BC/CB model
    zage(t,ttg)

    # Get number of redshift steps to use
    for kz in range(1,12000):
        z=zk(kz,zf)
        if z < zf:
            jl = kz
        else:
            jl=jl+1
            break

    # Get SED corresponding to redshift z = 0 (look-back-time = t0 = 0)
    t0, y0 = zsed(0.,t,f)

    # Compute magnitudes
    # print()
    # print(' Computing magnitudes/colors...')
    for kz in range(1,jl+1):
        z=zk(kz,zf)

        # Compute cosmological distance modulus.
        dm=cl.dismod(h,q,z)

        # Get SED corresponding to redshift z (look-back-time = tl)
        tl, yz = zsed(z,t,f)
        tz=ttg-tl

        # Read filter file and build filter arrays for redshift z = 0
        bt.ffd = fl.qfilters(w,kf,0.,0,1)

        # Compute magnitudes for y0 and yz at redshift 0
        f1_ne_0   = fl.VEGAmag(0,w,y0,0)
        ABf1_ne_0 =   fl.ABmag(0,w,y0,0)
        f1_rf     = fl.VEGAmag(0,w,yz,0)
        ABf1_rf   =   fl.ABmag(0,w,yz,0)
        f1_ev_0   = f1_rf
        if imag==2:
            f2_ne_0   = fl.VEGAmag(1,w,y0,0)
            ABf2_ne_0 =   fl.ABmag(1,w,y0,0)
            f2_rf     = fl.VEGAmag(1,w,yz,0)
            ABf2_rf   =   fl.ABmag(1,w,yz,0)
            f2_ev_0   = f2_rf

        # Build filter arrays for redshift z
        bt.ffd = fl.qfilters(w,kf,0.,z,1)

        # Compute absolute magnitudes in both filters for redshift z
        f1_ne     = fl.VEGAmag(0,w,y0,z)
        ABf1_ne   =   fl.ABmag(0,w,y0,z)
        f1_ev     = fl.VEGAmag(0,w,yz,z)
        ABf1_ev   =   fl.ABmag(0,w,yz,z)
        ABf1_ev_0 = ABf1_rf
        if imag==2:
            f2_ne     = fl.VEGAmag(1,w,y0,z)
            ABf2_ne   =   fl.ABmag(1,w,y0,z)
            f2_ev     = fl.VEGAmag(1,w,yz,z)
            ABf2_ev   =   fl.ABmag(1,w,yz,z)
            ABf2_ev_0 = ABf2_rf

        # Compute apparent magnitude
        f1mag    = f1_ev   + dm
        f1magAB  = ABf1_ev + dm
        # Compute k-correction
        #   For z=0 sed
        k1_corr_0   = f1_ne   - f1_ne_0
        k1_ABcorr_0 = ABf1_ne - ABf1_ne_0
        #   For evolved sed
        k1_corr_z   = f1_ev   - f1_rf
        k1_ABcorr_z = ABf1_ev - ABf1_rf
        # Compute (e+k)-correction
        ek1_corr   = f1_ev   - f1_ne_0
        ek1_ABcorr = ABf1_ev - ABf1_ne_0

        if imag==2:
            # Compute apparent magnitude
            f2mag    = f2_ev   + dm
            f2magAB  = ABf2_ev + dm
            # Compute k-correction
            #   For z=0 sed
            k2_corr_0   = f2_ne   - f2_ne_0
            k2_ABcorr_0 = ABf2_ne - ABf2_ne_0
            #   For evolved sed
            k2_corr_z   = f2_ev   - f2_rf
            k2_ABcorr_z = ABf2_ev - ABf2_rf
            # Compute (e+k)-correction
            ek2_corr   = f2_ev   - f2_ne_0
            ek2_ABcorr = ABf2_ev - ABf2_ne_0

            # Compute color kf(1) - kf(2) in Vega system
            c12_ev   = f1_ev   - f2_ev
            c12_rf   = f1_rf   - f2_rf
            c12_ne   = f1_ne   - f2_ne
            c12_ne_0 = f1_ne_0 - f2_ne_0
            c12_ev_0 = f1_ev_0 - f2_ev_0
            # Compute color kf(1) - kf(2) in AB system
            ABc12_ev   = ABf1_ev   - ABf2_ev
            ABc12_rf   = ABf1_rf   - ABf2_rf
            ABc12_ne   = ABf1_ne   - ABf2_ne
            ABc12_ne_0 = ABf1_ne_0 - ABf2_ne_0
            ABc12_ev_0 = ABf1_ev_0 - ABf2_ev_0
            # Compute color k-correction
            k12_corr_0   = c12_ne   - c12_ne_0
            k12_corr_z   = c12_ev   - c12_ev_0
            k12_ABcorr_0 = ABc12_ne - ABc12_ne_0
            k12_ABcorr_z = ABc12_ev - ABc12_ev_0
            # Compute color (e+k)-corrections
            ek12_corr   = c12_ev   - c12_ne_0
            ek12_ABcorr = ABc12_ev - ABc12_ne_0

        # Write results
        l=f'{z:5.3f}{tl:8.3f}{tz:8.3f}{dm:8.3f}   {f1_rf:9.4f}{f1_ne:9.4f}{f1_ev:9.4f}{f1mag:9.4f}   {ABf1_rf:9.4f}{ABf1_ne:9.4f}{ABf1_ev:9.4f}{f1magAB:9.4f}   {ek1_ABcorr:9.4f}{k1_ABcorr_z:9.4f}{k1_ABcorr_0:9.4f}\n'
        if z==0:
            l='#'+l
        else:
            l=' '+l
        o1.write(l)
        if imag==2:
            l=f'{z:5.3f}{tl:8.3f}{tz:8.3f}{dm:8.3f}   {f2_rf:9.4f}{f2_ne:9.4f}{f2_ev:9.4f}{f2mag:9.4f}   {ABf2_rf:9.4f}{ABf2_ne:9.4f}{ABf2_ev:9.4f}{f2magAB:9.4f}   {ek2_ABcorr:9.4f}{k2_ABcorr_z:9.4f}{k2_ABcorr_0:9.4f}\n'
            k=f'{z:5.3f}{tl:8.3f}{tz:8.3f}{dm:8.3f}   {c12_rf:9.4f}{c12_ne:9.4f}{c12_ev:9.4f}   {ABc12_rf:9.4f}{ABc12_ne:9.4f}{ABc12_ev:9.4f}   {ek12_ABcorr:9.4f}{k12_ABcorr_z:9.4f}{k12_ABcorr_0:9.4f}\n'
            if z==0:
                l='#'+l
            else:
                l=' '+l
            k=' '+k
            o2.write(l)
            o3.write(k)
        p = cn.percent(kz,jl,' Running code CM_EVOLUTION',p)
        if z >= zf:
            print()
            break

    # Report output files
    o1.close()
    print(h1)
    if imag==2:
        o2.close()
        o3.close()
        print(h2)
        print(h3)
    print()

def newp(file1):
    # Updates specific parameter
    global h, omega, omega_lambda, clambda, q, tu, ttg, zf
    print()
    print(' COSMOLOGY: Modify values of cosmological parameters if desired')
    print()
    print(' Enter redshift (z) and filter number (N) to compute magnitude evolution (to modify parameters in Set N, enter below z = -N (q to quit)')
    print()
    print('      Set 1 -> BC/CB model in file:      '   + bs.lnam(file1))                                 #print(' Set 1 -> BC/CB model in file:       ' + bs.lnam(file1))
    print('      Set 2 -> Hubble constant H0:       '   + str('{:.2f}'.format(h)) + ' Km/s/Mpc')          #print(' Set 2 -> Hubble constant H0:        70.00 Km/s/Mpc')
    print('               Omega:                     '  + str('{:.2f}'.format(omega)))                    #print('          Omega:                      0.30')
    print('               Omega_lambda:              '  + str('{:.2f}'.format(omega_lambda)))             #print('          Omega_lambda:               0.70')
    print('                 Lambda:                  '  + str('{:.2E}'.format(clambda)))                  #print('            Lambda:                   4.76E-05')
    print('                 Parameter q:            '   + str('{:.2f}'.format(q)))                        #print('            Parameter q:             -0.55')
    print('                 Age of the universe:    '   + str('{:.3f}'.format(tu)) + ' Gyr')              #print('            Age of the universe:     13.476 Gyr')
    print('      Set 3 -> Age of galaxy today:      '   + str('{:.3f}'.format(ttg)) + ' Gyr')             #print(' Set 3 -> Age of galaxy today:       13.000 Gyr')
    print('                 z of galaxy formation:  '   + str('{:.3f}'.format(zf)))                       #print('            z of galaxy formation:    9.842')
    print()
    z = input(' Enter z, filter number [ 0.000, 15] = ')
    return z

def cosmopar(file1):
    # Set and/or modify parameters if needed

    # Default values
    global h, omega, omega_lambda, clambda, q, tu, ttg, zf
    h, q, clambda, tu, ttg, zf, omega, omega_lambda = bt.cosmol

    while True:
       #cl.cosmol()
        z = newp(file1)
        if z == '-1':
           #s = updp(1,file1)
            file1 = updp(1,file1)
        elif z == '-2':
            s = updp(2,file1)
        elif z == '-3':
            s = updp(3,file1)
        elif z == 'q':
            print()
            return -1,15,file1
        else:
            if (z == ''):
                z = '0,15'
            s = z.replace(" ", ",")
            s = s.split(",")
            s = [x for x in s]
            z, n = float(s[0]), int(s[1])
            if z < 0.:
                z = 0.
            if n <= 0:
                n=15
            s = str('{:.4f}'.format(z)) + ',' + str('{:d}'.format(n))
            kf = [n]
            break
    return z,kf,file1

def magz(w,z,kf,y0,yz,tl,tz,dm,file1):
    global h, omega, omega_lambda, clambda, q, tu, ttg, zf

   # Compute absolute magnitudes in filter kf (Vega and AB systems)
    for i in range(len(kf)):

        # Compute absolute magnitudes in filter kf at z = 0

        # Read filter file and build filter arrays at z = 0
        bt.iread = True	# Read filter file
        bt.ffd   = fl.qfilters(w,kf,0.,0,1)

        # Vega magnitude of evolved sed at redshift 0
        f_rf = fl.VEGAmag(i,w,yz,0.)

        # Vega magnitude of z = 0 sed at redshift 0
        f_ne_0 = fl.VEGAmag(i,w,y0,0.)

        # AB magnitude of evolved sed at redshift 0
        ABf_rf = fl.ABmag(i,w,yz,0)

        # AB magnitude of z = 0 sed at redshift 0
        ABf_ne_0 = fl.ABmag(i,w,y0,0)

        # Compute absolute magnitudes in filter kf at z

        # Build filter arrays at z = 0
        bt.ffd = fl.qfilters(w,kf,0.,z,1)

        # Vega magnitude of z = 0 sed at redshift z
        f_ne = fl.VEGAmag(i,w,y0,z)

        # Vega magnitude of evolved sed at redshift z
        f_ev = fl.VEGAmag(i,w,yz,z)

        # AB magnitude of z = 0 sed at redshift z
        ABf_ne = fl.ABmag(i,w,y0,z)

        # AB magnitude of evolved sed at redshift z
        ABf_ev = fl.ABmag(i,w,yz,z)

        # Compute apparent magnitudes
        fmag_NG = f_ev   + dm
        fmagAB  = ABf_ev + dm

        # Compute k-corrections:
        #   For z=0 sed
        k_corr_0   = f_ne   - f_ne_0
        k_ABcorr_0 = ABf_ne - ABf_ne_0
        #   For evolved sed
        k_corr_z   = f_ev   - f_rf
        k_ABcorr_z = ABf_ev - ABf_rf

        # Compute (e+k)-corrections
        ek_corr   = f_ev   - f_ne_0
        ek_ABcorr = ABf_ev - ABf_ne_0
 
        # Print results
        n = kf[i]
        h, q, clambda, tu, ttg, zf, omega, omega_lambda = bt.cosmol2
        print()
        print('      Set 1 -> BC/CB model in file:      '   + bs.lnam(file1))                                 #print(' Set 1 -> BC/CB model in file:       ' + bs.lnam(file1))
        print('      Set 2 -> Hubble constant H0:       '   + str('{:.2f}'.format(h)) + ' Km/s/Mpc')          #print(' Set 2 -> Hubble constant H0:        70.00 Km/s/Mpc')
        print('               Omega:                     '  + str('{:.2f}'.format(omega)))                    #print('          Omega:                      0.30')
        print('               Omega_lambda:              '  + str('{:.2f}'.format(omega_lambda)))             #print('          Omega_lambda:               0.70')
        print('                 Lambda:                  '  + str('{:.2E}'.format(clambda)))                  #print('            Lambda:                   4.76E-05')
        print('                 Parameter q:            '   + str('{:.2f}'.format(q)))                        #print('            Parameter q:             -0.55')
        print('                 Age of the universe:    '   + str('{:.3f}'.format(tu)) + ' Gyr')              #print('            Age of the universe:     13.476 Gyr')
        print('      Set 3 -> Age of galaxy today:      '   + str('{:.3f}'.format(ttg)) + ' Gyr')             #print(' Set 3 -> Age of galaxy today:       13.000 Gyr')
        print('                 z of galaxy formation:  '   + str('{:.3f}'.format(zf)))                       #print('            z of galaxy formation:    9.842')
        print()
        print(' GALAXEV model sed in file: ' + bs.lnam(file1))
        print()
        print(' Meaning of printed variables (see Tables 5 and 6 of bc03.pdf for details):')
        print('   z       = redshift')
        print('   ltt     = light travel time from z=0 to z in Gyr')
        print('   tz      = age of galaxy at redshit = z')
        print('   dm      = cosmological distance modulus in magnitude units')
        print('   M_rf    = M_rf(z) in filter kf')
        print('   M_ne    = M_ne(z) in filter kf')
        print('   M_ev    = M_ev(z) in filter kf')
        print('   m_ev    = apparent magnitude M_ev(z) + dm(z) in filter kf')
        print('   e+kcor  = (e+k)-correction in filter kf')
        print('   k_cor_z = k-correction in filter kf computed from sed at redshift z')
        print('   k_cor_0 = k-correction in filter kf computed from sed at redshift 0')
        print('   Z.P.    = Zeropoint definition, either Vega or AB system')
        print()
        print(' For filter No. ' + str(n) + ' = ' + fl.filterid(n))
        print('    z   ltt(Gyr) tz(Gyr)   dm      M_rf     M_ne     M_ev     m_ev    e+kcor   k_cor_z  k_cor_0  P.S.')
        l = f'{z:7.3f}{tl:8.3f}{tz:8.3f}{dm:8.3f}{f_rf:9.4f}{f_ne:9.4f}{f_ev:9.4f}{fmag_NG:9.4f}{ek_corr:9.4f}{k_corr_z:9.4f}{k_corr_0:9.4f}' + '   Vega'
        print(l)
        l = ''
        for i in range(31): l=l+' '
        l = l + f'{ABf_rf:9.4f}{ABf_ne:9.4f}{ABf_ev:9.4f}{fmagAB:9.4f}{ek_ABcorr:9.4f}{k_ABcorr_z:9.4f}{k_ABcorr_0:9.4f}' + '    AB'
        print(l)
        l = ' '
        for i in range(62): l=l+'-'
        l=l+'^'
        for i in range(38): l=l+'-'
        print(l)
 
def zmag(file1):
    # python version of zmag code

    global f, t, fd, iread

    # Init variables
    bt.iread = True     # Read filter file only on first call
    p=-1                # Needed by 'percent' function on first call

    # Ask for input parameters to run zmag program
    bs.bchead()
    f1 = input(' BC_GALAXEV SSP sed in file [' + bs.lnam(file1) + '] = ')
    if len(f1) <= 0:
        file1 = bs.lnam(file1)
    else:
        file1 = bs.zrep(file1,f1)

    # Read BC/CB model in fits file
    w,f,t,e,*nouse = bc.bcfits(file1)
    hd.multhead(0,file1,t,w,e)

    while True:
        # Get parameters
        z, kf, file2 = cosmopar(file1)
        if z < 0:
            break
        bt.cosmol2 = h, q, clambda, tu, ttg, zf, omega, omega_lambda

        # Read BC/CB model in fits file
        if file2 != file1:
            w,f,t,e,*nouse = bc.bcfits(file2)
            hd.multhead(0,file2,t,w,e)
            file1 = file2

        # Get z corresponding to each look-back-time in BC/CB model
        zage(t,ttg)

        # Get SED corresponding to redshift z = 0 (look-back-time = t0 = 0)
        t0, y0 = zsed(0.,t,f)

        # Get SED corresponding to redshift z (look-back-time = tl)
        tl, yz = zsed(z,t,f)

        # Galaxy age at look-back-time tl
        tz = ttg-tl

        # Compute cosmological distance modulus
        dm=cl.dismod(h,q,z)

        # Compute magnitudes vs. z
        magz(w,z,kf,y0,yz,tl,tz,dm,file1)

def filesys(inp):
    # File name (asssumed fits file)
    inp[0] = bs.lnam(inp[0])
    file = inp[0]
    if '.fits' not in file:
        file = file + '.fits'
    # Check for user defined filters or predefined filter set
    if 'user' in inp[1][0]:
        if 'VEGA' in inp[1][1]:
            fds = 'user.VEGAmag'
        elif 'AB' in inp[1][1]:
            fds = 'user.ABmag'
        elif 'ST' in inp[1][1]:
            fds = 'user.STmag'
        i = [eval(i) for i in inp[1][2:]]
        bt.usr = i
        inp[1] = fds
    elif 'jplus' in inp or 'JPLUS' in inp:
        fds = 'jplus'
        i   = np.arange(296, 308, dtype=int)
    elif 'jpas' in inp or 'JPAS' in inp:
        fds = 'jpas'
        i   = np.arange(312, 372, dtype=int)
    elif 'sdss' in inp or 'SDSS' in inp:
        fds = 'sdss'
        i   = np.arange(120, 125, dtype=int)
    elif 'stmag' in inp or 'STmag' in inp:
        fds = 'stmag'
        i   = np.concatenate( ([308, 309, 219, 220, 310, 222, 311, 224, 211], np.arange(225, 230, dtype=int), np.arange(276, 293, dtype=int)) )
    elif 'johnson' in inp or 'JOHNSON' in inp or 'john' in inp or 'JOHN' in inp:
        fds = 'johnson'
        i   = [12, 13, 14, 15, 84, 85, 32, 33, 34, 35, 36, 55, 56, 57, 88]
    return file,fds,i

def getp(a,inp):
    # Assigns parameters read from the command line
    # Finds character 'a' in inp
    if a=='-o':
        b=True
        return b
    # Avoid arguments after '-o'
    c=[]
    for i in range(2,len(inp)):
        if inp[i] == '-o':
            break
        c.append(inp[i])
    # Search in string c
    b=-1.
    for i in range(len(c)):
        if a in c[i]:
            b = c[i]
            b = b.replace(a,'')
            b = b.replace('=','')
            b = float(b)
    b = max(0.,b)
    return b

def rf_phot():
    # Computes magnitudes for a galaxy observed at redshift z at any age
    global h, omega, omega_lambda, clambda, q, tu, ttg, zf

    # Check number of arguments
    inp = bt.rf
    if len(inp) < 2:
        print('Sample usage: rf_phot.py file.fits filters [z=0.1] [av=0.5] [t=10] [-o] [output_file_name]')
        sys.exit()

    # Check if fits is in file name
#   if '.fits' not in inp[0]:
#       inp[0]  = inp[0] + '.fits'

    # Interpret command line parameters
    file,fds,i = filesys(inp)       # Get file name and filter system and filters to use
    tx = getp('t',inp)              # Check for 't' in argument list
    av = getp('av',inp)             # Check for 'av' in argument list
    z  = getp('z',inp)              # Check for 'z' in argument list

    # Compute cosmological distance modulus.
    h, q, clambda, tu, ttg, zf, omega, omega_lambda = cl.cosmol()
    dm=cl.dismod(h,q,z)

    # Read BC/CB model in fits file
    w,f,t,*nouse = bc.bcfits(file)

    # Open output files
    of = hd.opfiles(inp,z,av)

    # Build time scale
    if tx <= 0:
        # Do all time steps
        # td = t
        # Do following time steps (in Gyr)
        td = [ 0.0005,  0.0008,  0.0010,  0.0012,  0.0014,  0.0016,  0.0018,  0.0020,  0.0025,  0.0030, 0.0035, 0.0040, 0.0045, 0.0050, 0.0060, \
               0.0070,  0.0080,  0.0090,  0.0100,  0.0120,  0.0140,  0.0160,  0.0180,  0.0200,  0.0220, 0.0240, 0.0260, 0.0280, 0.0300, 0.0320, \
               0.0340,  0.0360,  0.0380,  0.0400,  0.0425,  0.0450,  0.0475,  0.0500,  0.0535,  0.0570, 0.0600, 0.0700, 0.0800, 0.0900, 0.1000, \
               0.1200,  0.1400,  0.1600,  0.1800,  0.2000,  0.2250,  0.2500,  0.2750,  0.3000,  0.3500, 0.4000, 0.4500, 0.5000, 0.5500, 0.6000, \
               0.6500,  0.7000,  0.7500,  0.8000,  0.8500,  0.9000,  1.0000,  1.2000,  1.4000,  1.6000, 1.8000, 2.0000, 2.2000, 2.4000, 2.5000, \
               2.6000,  2.8000,  3.0000,  3.5000,  4.0000,  4.5000,  5.0000,  5.5000,  6.0000,  6.5000, 7.0000, 7.5000, 8.0000, 8.5000, 9.0000, \
               9.5000, 10.0000, 10.5000, 11.0000, 11.5000, 12.0000, 12.5000, 13.0000, 13.5000, 14.0000 ]
    else:
        # Do time step tx
        td=[]
        td.append(tx)

    print()
    print(bs.rp(' Running rf_phot code (may take a while)'))

    # Read filters and add model wavelength points at corresponding z
    bt.iread = True     # Read filter file on first call
    p=-1                # Needed by 'percent' function on first call

    # Build filter arrays for redshift z
    bt.ffd   = fl.qfilters(w,i,av,z,1)

    # Compute and report requested magnitudes
    jl = len(td)
    for k in range(jl):
        # Get sed corresponding to age a
        a = td[k]
        y = bc.ft(a,f,t)
        # Compute magnitudes
        if fds == 'sdss' or fds == 'jplus' or fds == 'jpas' or fds == 'user.ABmag':
            # Compute ABmag and Flux in Jansky
            fl.ABphot(w,y,z,dm,0.,0.,av,a,of,0)
        elif fds == 'stmag' or fds == 'user.STmag':
            # Compute STmag and Flux in STflux
            fl.STphot(w,y,z,dm,0.,0.,av,a,of,0)
        elif fds == 'johnson' or fds =='user.VEGAmag':
            fl.VEGAphot(w,y,z,dm,0.,0.,av,a,of,0)
        else:
            print('Unknown filter system:',fds)
            quit()
        p = cn.percent(k,jl,' Running code RF_PHOT',p)

    print()
    if of:
        print()
        print(' Output written to file(s):',of[3].replace('//','/'))
        print('                           ',of[4].replace('//','/'))

def of_phot():
    # Computes observer frame magnitudes vs z, Av
    global h, omega, omega_lambda, clambda, q, tu, ttg, zf

    # Check number of arguments
    inp = bt.of
    if len(inp) < 2:
        print('Sample usage: of_phot.py file.fits filters [av=0.5]')
        sys.exit()

    # Check if fits is in file name
    if '.fits' not in inp[0]:
        inp[0]  = inp[0] + '.fits'

    # Interpret command line parameters
    file,fds,i = filesys(inp)       # Get file name and filter system and filters to use
    av = getp('av',inp)             # Check for 'av' in argument list

    # Get values of cosmological parameters
    if not jupy:
        # Default
        h, q, clambda, tu, ttg, zf, omega, omega_lambda = bt.cosmol
        bt.cosmol2 = h, q, clambda, tu, ttg, zf, omega, omega_lambda
    else:
        # May have been changed by user
        h, q, clambda, tu, ttg, zf, omega, omega_lambda = bt.cosmol2

    print()
    print(bs.rp(' Running of_phot code (may take a while)'))

    # Read BC/CB model in fits file
    w,f,t,*nouse = bc.bcfits(file)

    # Open output files
    of = hd.oqfiles(inp,av)

    # Get z corresponding to each look-back-time in BC/CB model
    zage(t,ttg)

    # Init variables
    bt.iread = True     # Read filter file only on first call
    p=-1                # Needed by 'percent' function on first call

    # Get number of redshift steps to use
    for kz in range(1,12000):
        z=zk(kz,zf)
        if z < zf:
            jl = kz
        else:
            jl=jl+1
            break

    # Compute magnitudes
    for kz in range(1,jl+1):
        z=zk(kz,zf)

        # Compute cosmological distance modulus.
        dm=cl.dismod(h,q,z)

        # Get SED corresponding to redshift z (look-back-time = tl)
        tl, yz = zsed(z,t,f)
        tz=ttg-tl

        # Build filter arrays for redshift z
        bt.ffd = fl.qfilters(w,i,av,z,1)

        # Compute magnitudes
        if fds == 'sdss' or fds == 'jplus' or fds == 'jpas' or fds == 'user.ABmag':
            # Compute ABmag and Flux in Jansky
            fl.ABphot(w,yz,z,dm,tl,tz,av,z,of,1)
        elif fds == 'stmag' or fds == 'user.STmag':
            # Compute ABmag and Flux in Jansky
            fl.STphot(w,yz,z,dm,tl,tz,av,z,of,1)
        elif fds == 'johnson' or fds == 'user.VEGAmag':
            fl.VEGAphot(w,yz,z,dm,tl,tz,av,z,of,1)
        else:
            print('Unknown filter system:',fds)

        p = cn.percent(kz,jl,' Running code OF_PHOT',p)

        if z >= zf:
            break
    print()
    print(' Output written to file(s):',of[3].replace('//','/'))
    print()

def krfp(file):

    # Check if parameters have been entered in command line
    if isinstance(file, str):
        file = [file, 'user']
    elif isinstance(file, list) and len(file) > 0 and 'user' not in file:
        bt.rf = file
        return

    # Ask for input parameters to run rf_phot program
    bs.bchead()
    f1 = input(' BC_GALAXEV SSP sed in file [' + bs.lnam(file) + '] = ')
    if len(f1) <= 0:
        file = bs.lnam(file)
    else:
        file = bs.zrep(file,f1)
    if '.fits' not in file:
        file = file + '.fits'

    # Select filter sets
    print()
    print(' The following filter sets are defined:')
    print('     Set: 1 - SDSS (ugriz filters)')
    print('          2 - JPAS (60 J-PAS filters)')
    print('          3 - JPLUS (12 J-PLUS filters)')
    print('          4 - STmag (14 HST-ACS, 8 JWST-NIRCAM, 9 JWST-MIRI filters)')
    print('          5 - Johnson (UBVRIJKL RcIc PalJHK K'' filters)')
    print('          6 - enter filters of your choice')
    print('          q - quit')
   #print()
    i = input(' choice = ')
    if i == '2':
        fds = 'jpas'
        mag = '.ABmag'
    elif i == '3':
        fds = 'jplus'
        mag = '.ABmag'
    elif i == '4':
        fds = 'stmag'
        mag = '.STmag'
    elif i == '5':
        fds = 'johnson'
        mag = '.VEGAmag'
    elif i == '6':
        gds = input(' Enter number of selected filters = ')
        mag = input(' Desired magnitude: 0:VEGA, 1:AB, 2:ST = ')
        fds = 'user'
        if mag == '0':
            mag = '.VEGAmag'
        elif mag == '2':
            mag = '.STmag'
        else:
            mag = '.ABmag'
        # User selected filters
        if len(gds) > 0:
            gds = gds.replace(',',' ')
            gds = 'user ' + mag + ' ' + gds
            gds = gds.split()
        else:
            sys.exit()
    elif i == 'q':
        sys.exit()
    else:
        fds = 'sdss'
        mag = '.ABmag'

    # Ask for redshift to compute magnitudes
   #print()
    zs = input(' Compute rest frame magnitudes at redshift z [0] = ')
    if zs == '':
        z = 0.
        b = 'z0.0'
    else:
        z = float(zs)
        b = 'z' + zs

    # Ask for extinction Av
   #print()
    es = input(' Redden SEDs assuming Av [0] = ')
    if es == '':
        av = 0.
        a  = 'av0.0'
    else:
        av = float(es)
        a  = 'av' + es

    # Build output file name
    n = file.replace('fits','rf_phot.' + fds)
    n = n + ('.z' + f'{z:6.3f}' + '.Av' + f'{av:6.3f}').replace(' ','')
    n = n + mag
   #print()
    of = input(' Output file name [' + n + '] = ')

    # Build command line arguments to run code
    if fds == 'user':
        fds = gds
    if of == '':
        bt.rf = [file, fds, b, a, '-o']
    else:
        bt.rf = [file, fds, b, a, '-o', of]

def kofp(file):

    # Check if parameters have been entered in command line
    if isinstance(file, str):
        file = [file, 'user']
    elif isinstance(file, list) and len(file) > 0 and 'user' not in file:
        bt.of = file
        return

    # Ask for input parameters to run of_phot program
    bs.bchead()
    f1 = input(' BC_GALAXEV SSP sed in file [' + bs.lnam(file) + '] = ')
    if len(f1) <= 0:
        file = bs.lnam(file)
    else:
        file = bs.zrep(file,f1)

    # Select filter sets
    print()
    print(' The following filter sets are defined:')
    print('     Set: 1 - SDSS (ugriz filters)')
    print('          2 - JPAS (60 J-PAS filters)')
    print('          3 - JPLUS (12 J-PLUS filters)')
    print('          4 - STmag (14 HST-ACS, 8 JWST-NIRCAM, 9 JWST-MIRI filters)')
    print('          5 - Johnson (UBVRIJKL RcIc PalJHK K'' filters)')
    print('          6 - enter filters of your choice')
    print('          q - quit')
   #print()
    i = input(' choice = ')
    if i == '2':
        fds = 'jpas'
        mag = '.ABmag'
    elif i == '3':
        fds = 'jplus'
        mag = '.ABmag'
    elif i == '4':
        fds = 'stmag'
        mag = '.STmag'
    elif i == '5':
        fds = 'johnson'
        mag = '.VEGAmag'
    elif i == '6':
        gds = input(' Enter number of selected filters = ')
        mag = input(' Desired magnitude: 0:VEGA, 1:AB, 2:ST = ')
        fds = 'user'
        if mag == '0':
            mag = '.VEGAmag'
        elif mag == '2':
            mag = '.STmag'
        else:
            mag = '.ABmag'
        # User selected filters
        if len(gds) > 0:
            gds = gds.replace(',',' ')
            gds = 'user ' + mag + ' ' + gds
            gds = gds.split()
        else:
            sys.exit()
    elif i == 'q':
        sys.exit()
    else:
        fds = 'sdss'
        mag = '.ABmag'

    # Ask for extinction Av
   #print()
    es = input(' Redden SEDs assuming Av [0] = ')
    if es == '':
        av = 0.
        a  = 'av0.0'
    else:
        av = float(es)
        a  = 'av' + es

    # Build output file name
    n = file.replace('fits','of_phot.' + fds)
    n = n + ('.Av' + f'{av:6.3f}').replace(' ','')
    n = n + mag
   #print()
    of = input(' Output file name [' + n + '] = ')

    # Build command line arguments to run code
    if fds == 'user':
        fds = gds
    if of == '':
        bt.of = [file, fds, a, '-o']
    else:
        bt.of = [file, fds, a, '-o', of]

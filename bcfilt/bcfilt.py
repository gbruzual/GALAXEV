import os
import time
import warnings
import numpy as np
import builtins as bt
from   astropy.io import fits
from   astropy.utils.exceptions import AstropyWarning
from   astropy.table import Table
import basics.basics as bs
import bchead.bchead as hd
import bcfits.bcfits as bf

def qfilters(w,g,av,z,o):
    global p, r
    # Read filters
    if bt.iread:
        p,r = readfilters()
        bt.iread = False

    # If o == 0, add filter wavelength points to SED wavelength points (as in fortran code)
    #            this option is preferred for low resolution spectra

    # If o != 0, do not add filter wavelength points to SED wavelength
    #            this option is faster in python and works well for high resolution spectra

    if o == 0:
        # Prepare filters for usage
        fw = []
        fo = []
        fr = []
        fn = []
        ip = []
        i1 = []
        i2 = []
        j1 = []
        j2 = []
        zp = []
        m = 0
        for j in range (len(g)):
            # Add filter wavelength points to SED wavelength
            k, mat, ro, l1, l2, zo = ofixfilters(abs(g[j]),p,r,w,z)
            ip.append(k)
            i1.append(m)
            m = m + k
            i2.append(m-1)
            j1.append(l1)
            j2.append(l2)
            zp.append(zo)
            for l in range(k):
                fn.append(g[j])
                fw.append(mat[l][0])
                fr.append(mat[l][1])
            fo = fr
        return ip, i1, i2, j1, j2, fw, fr, fo, fn, zp

    elif o != 0:
        # Read reddening law
        xd, yd = q_ccm(av)

        # Prepare filters for usage
        fw = []
        fo = []
        fr = []
        fn = []
        ip = []
        i1 = []
        i2 = []
        j1 = []
        j2 = []
        zp = []
        m = 0
       #start = time.time()
        for j in range (len(g)):
            # Do not add filter wavelength points to SED wavelength
            k, xt, ro, rt, l1, l2, zo = qfixfilters(abs(g[j]),p,r,w,z,av,xd,yd)
            ip.append(k)
            i1.append(m)
            m = m + k
            i2.append(m-1)
            j1.append(l1)
            j2.append(l2)
            zp.append(zo)
            gj=[g[j]]*k
            fn.extend(gj)
            fw.extend(xt)
            fo.extend(ro)
            fr.extend(rt)
           #for l in range(k):
           #    fn.append(g[j])
           #    fw.append(xt[l])
           #    fo.append(ro[l])
           #    fr.append(rt[l])
       #end = time.time()
       #print(end-start)
        return ip, i1, i2, j1, j2, fw, fr, fo, fn, zp

def qfixfilters(i,p,r,w,z,av,xd,yd):
    # >>> use filter only at SED wavelength points <<<

    # Read filter I
    xf,rf,nf,zo = getfilter(i,p,r)

    # Shift filter wavelength to the blue by (1+z) as in function FILTER_N
    z1 = 1.+z
    for i in range(len(xf)):
        xf[i] = xf[i]/z1

    # Interpolate blue-shifted filter transmission function at the SED wavelengths
    # Store original transmission function in array ro
    mp = len(xf)-1
    j1 = np.searchsorted(w, xf[0])
    j2 = np.searchsorted(w, xf[mp])
    xt = []
    rt = []
    ro = []
    for i in range (j1,j2):
            y = np.interp(w[i], xf, rf, left=0., right=0.)
           #y = linear(w[i], xf, rf, left=0., right=0.)
            xt.append(w[i])
            ro.append(y)
            # To save computer time, apply reddening correction to the filter transmission
            if av > 0:
                yi = np.interp(w[i], xd, yd, left=yd[0], right=yd[len(yd)-1])
                y  = 10.**(-0.4*av*yi)*y
            rt.append(y)
    return len(xt), xt, ro, rt, j1, j2, zo

def ofixfilters(i,p,r,w,z):
    # >>> Follow same procedure as in fortran code (function FILTER_N) <<<

    # From FILTER_N:
    #   Improved version that samples correctly the s.e.d. and the
    #   filter response function at all points available in these arrays.

    #   The current version of the F_MEAN routine requires the filter
    #   response function to be interpolated at each point in the
    #   sed which is not a point in the filter response. All points in
    #   the filter response are also used. The format of the binary
    #   file was changed to reduce its size. Filters are stored
    #   sequentially in r(i,k), k=1 => wavelength, k=2 => response
    #   function. Information about starting point and number of
    #   points per filter is kept in file, as well as filter id label.

    # Read filter I and shift filter wavelength to the blue by (1+z)
    #   Arrays (xt,rt) are required for np.interp to interpolate correctly,
    #   it does not work correctly using arrays (xf,rf) which are modified
    #   by the append command after each interpolation.
    xf,rf,nf,zp = getfilter(i,p,r)					# ! In FILTER_N:
    z1 = 1.+z
    xt = []                    					# ! Shift filter wavelength by (1+z)
    rt = []                    					# m=0
    ro = []                    					# m=0
    for i in range(len(xf)):   					# do k=ni(i),nl(i)
        xf[i] = xf[i]/z1       					# m=m+1
        xt.append(xf[i])       					# xf(m)=r(k,1)/z1
        rt.append(rf[i])       					# rf(m)=r(k,2)
        ro.append(rf[i])       					# rf(m)=r(k,2)
      								# enddo
    # Add wavelength points in the model sed (x,y) in the galaxy
    # restframe to array xf interpolating transmission function
    mp = len(xf)-1
    j1 = np.searchsorted(w, xf[0])
    j2 = np.searchsorted(w, xf[mp])				# ! In FILTER_N:
    for i in range (len(w)):					# ! Add wavelength points in the sed (x,y) in the galaxy restframe
        if w[i] > xf[mp]:                                       # l=0
            break                                               # do k=1,n
        elif w[i] >= xf[0]:                                     # ! xz=x(k)/z1    ! This statement is wrong. x(i) is already in the galaxy rest frame
            y = np.interp(w[i], xt, rt, left=0., right=0.)      # xz=x(k)         ! This is correct. Emma C. Lake noticed this error in March 2017. !!!!!!!!!!!
           #y = linear(w[i], xt, rt, left=0., right=0.
            xf.append(w[i])                                     # if (xz >= xf(np(i))) then
            rf.append(y)                                        #	goto 1
    # Build matrix to sort (xf,rf) pairs according to xf        # elseif (xz >= xf(1)) then
    rows = len(xf)                                              #	m=m+1
    cols = 2                                                    #	xf(m)=xz
    matrix = []                                                 #	! Interpolate shifted filter at shifted wavelength of sed
    for i in range(rows):                                       #	rf(m)=LINEAR(xz,xf,rf,np(i),l)
        row = []                                                # endif
        for j in range(cols):                                   # enddo
            if j==0:                                            # ! Sort (xf,rf) arrays according to xf
                row.append(xf[i])                               # 1 call sort2(m,xf,rf)
            else:
                row.append(rf[i])
        matrix.append(row)
    # Sort along the first axis
    matrix.sort()
    return rows, matrix, ro, j1, j2, zp

def readfilters():
    # Read filter file in fits table format
   #namf = os.environ.get('FILTERS') + '/FILTERRES.fits'
    namf = os.environ.get('FILTERF')
    namf = namf.replace('//','/')
    print()
    print(' Preparing filters in file: ' + namf)
    with fits.open(namf) as hdul:
        # Use Table to read the fits table
        n = Table.read(hdul,hdu=1)      # hdu = 1 => Filter record information
        r = Table.read(hdul,hdu=2)      # hdu = 2 => Filter transmission
    return n,r

def getfilter(i,p,r):
    # Returns response function for filter I
    n1 = p[i-1][1]
    n2 = p[i-1][2]
    zp = p[i-1][3]
    # Build array with response function
    w = []
    t = []
    for k in range(n1-1,n2):
        w.append(r[k][0])
        t.append(r[k][1])
    return w,t,len(w),zp

def q_ccm(av):
    # Reads extinction law from file ExtLaws.out
    if av > 0:
        # Read ExtLaw from ascii file
        # d = bf.nfile(os.environ.get('glxaux') + '/ExtLaws.out',1)
        # xd = []
        # yd = []
        # for i in range(len(d)):
        #   xd.append(d[i][0])
        #   yd.append(d[i][1])

        # Write fits file with q_ccm Extinction Law
        # namo = os.environ.get('glxaux') + '/ExtLaws.tmp'
        # namf = os.environ.get('glxaux') + '/ExtLaws.fits'
        # o = open(namo,'w')
        # o.write('#lambda       q_CCM\n')
        # for i in range (len(xd)):
        #  s = '%7.1f' % xd[i] + '%12.5f' % yd[i] + '\n'
        #  o.write(s)
        # o.close()
        # Create fits file
        # t = Table.read(namo, format='ascii')
        # d.write(namf, overwrite=True)
        # quit()

        # Read ExtLaw from fits file
        # namf = os.environ.get('glxaux') + '/ExtLaws.fits'
        # with fits.open(namf) as hdul:
        #     # Use Table to read the fits table
        #     n = Table.read(hdul,hdu=1)
        # xd = n['col1']
        # yd = n['col2']
        namf = os.environ.get('glxaux') + '/ExtLaw1.fits'
        with fits.open(namf) as hdul:
            # Use Table to read the fits table
            n = Table.read(hdul,hdu=1)
        xd = n['lambda']
        yd = n['q_CCM']
    else:
        xd=0
        yd=0
    return xd,yd

def qfilter_n(k,x,y,z):
    fd = bt.ffd
    i1 = fd[1] ; i2 = fd[2] ; j1 = fd[3] ; xlam = fd[5] ; rlam = fd[6]
    # ! In FILTER_N:
    # ! Compute filter_n = number of photons through filter
    # ! Compute number of photons, instead of total flux (Fukugita et al. 1996)
    # !   xlam() = filter wavelength
    # !   xf()   = filter wavelength shifted by (1.+z)
    # !   x()    = bc model wavelength
    # do j=pos_i(k),pos_f(k)
    # m=m+1
    # xf(m)=xlam(j)
    # rf(m)=z1*xlam(j)*rlam(j)*LINEAR(xf(m),x,y,n,l)
    # enddo
    # ! Factor (1+z) already included in lambda scale.
    # filtr_n=TRAPZ1(xf,rf,m)
    z1  = 1.+z
    xf  = []
    rfa = []
    m = 0
    for j in range(i1[k],i2[k]):
        m=m+1
        xa = xlam[j]   			# filter wavelength in the galaxy frame
        xf.append(xa)
        rfa.append(xa*z1*rlam[j]*y[j1[k]+m-1])
    filter_n = np.trapz(rfa,xf)
    return filter_n

def ofilter_n(k,x,y,z):
    fd = bt.ffd
    i1 = fd[1] ; i2 = fd[2] ; xlam = fd[5] ; rlam = fd[6]
    # ! In FILTER_N:
    # ! Compute filter_n = number of photons through filter
    # ! Compute number of photons, instead of total flux (Fukugita et al. 1996)
    # !   xlam() = filter wavelength
    # !   xf()   = filter wavelength shifted by (1.+z)
    # !   x()    = bc model wavelength
    # do j=pos_i(k),pos_f(k)
    # m=m+1
    # xf(m)=xlam(j)
    # rf(m)=z1*xlam(j)*rlam(j)*LINEAR(xf(m),x,y,n,l)
    # enddo
    # ! Factor (1+z) already included in lambda scale.
    # filtr_n=TRAPZ1(xf,rf,m)
    z1  = 1.+z
    xf  = []
    rfa = []
    for j in range(i1[k],i2[k]):
        xa = xlam[j]			# filter wavelength in the galaxy frame
        xf.append(xa)
        ya = np.interp(xlam[j], x, y, left=0., right=0.)
        rfa.append(z1*xa*rlam[j]*ya)
    filter_n = np.trapz(rfa,xf)
    return filter_n

def qf_mean(k,x,y,z):
    fd = bt.ffd
    i1 = fd[1] ; i2 = fd[2] ; j1 = fd[3] ; xlam = fd[5] ; rlam = fd[6] ; olam = fd[7]
    # ! In FILTER_N:
    # ! Compute f_mean = mean flux Fnu through filter (as required by AB mags)
    # ! Returns mean flux >>> (Fnu) <<< in ith filter (/Hz)
    # ! As: f_nu=int(dnu Fnu Rnu/h*nu)/int(dnu Rnu/h*nu)
    # ! ie: f_nu=int(dlm Flm Rlm lm / c)/int(dlm Rlm/lm)
    # !   xlam() = filter wavelength
    # !   xf()   = filter wavelength shifted by (1.+z)
    # !   x()    = bc model wavelength
    #
    # do j=pos_i(k),pos_f(k)
    # m=m+1
    # xf(m) =xlam(j)*z1       ! filter wavelength in detector's frame
    # rfa(m)=rlam(j)*xf(m)*LINEAR(xlam(j),x,y,n,l)
    # rfb(m)=rlam(j)/xf(m)
    # enddo
    #
    # f_mean=TRAPZ1(xf,rfa,m)/TRAPZ1(xf,rfb,m)
    # f_mean=f_mean/2.997925e+18              ! speed of light in A/sec
    # ! F(lambda)*dlambda = F[lambda/(1+z)]*dlambda/(1+z)
    # f_mean=f_mean/z1
    z1  = 1.+z
    xf  = []
    rfa = []
    rfb = []
    m = 0
    for j in range(i1[k],i2[k]):
        m=m+1
        xa = xlam[j]*z1			# filter wavelength in detector's frame
        xf.append(xa)
        rfa.append(rlam[j]*xa*y[j1[k]+m-1])
        rfb.append(olam[j]/xa)
    f_mean = np.trapz(rfa,xf) / np.trapz(rfb,xf)
    f_mean=f_mean/2.997925e+18		# speed of light in A/sec
    # F(lambda)*dlambda = F[lambda/(1+z)]*dlambda/(1+z)
    f_mean=f_mean/z1
    return f_mean

def qf_meanl(k,x,y,z):
    fd = bt.ffd
    i1 = fd[1] ; i2 = fd[2] ; j1 = fd[3] ; xlam = fd[5] ; rlam = fd[6] ; olam = fd[7]
    # ! In FILTER_N:
    # ! Compute f_meanl = mean flux Flambda through filter (as required by ST mags)
    # ! Returns mean flux >>> (Flambda) <<< in ith filter (/A)
    # ! As: f_lm=int(dlm Flm Rlm lm/h*c)/int(dlm Rlm lm/h*c)
    # ! ie: f_lm=int(dlm Flm Rlm lm)/int(dlm Rlm lm)
    # !   xlam() = filter wavelength
    # !   xf()   = filter wavelength shifted by (1.+z)
    # !   x()    = bc model wavelength
    #
    # do j=pos_i(k),pos_f(k)
    # m=m+1
    # xf(m) =xlam(j)*z1       ! wavelength in detector''s frame
    # rfa(m)=rlam(j)*xf(m)*LINEAR(xlam(j),x,y,n,l)
    # rfb(m)=rlam(j)*xf(m)
    # enddo
    #
    # f_meanl=TRAPZ1(xf,rfa,m)/TRAPZ1(xf,rfb,m)
    # f_meanl=f_meanl/z1
    z1  = 1.+z
    xf  = []
    rfa = []
    rfb = []
    m = 0
    for j in range(i1[k],i2[k]):
        m=m+1
        xa = xlam[j]*z1			# filter wavelength in detector's frame
        xf.append(xa)
        rfa.append(rlam[j]*xa*y[j1[k]+m-1])
        rfb.append(olam[j]*xa)
    f_meanl = np.trapz(rfa,xf) / np.trapz(rfb,xf)
    f_meanl=f_meanl/z1
    return f_meanl

def absmag(x,y,kf):
    # Compute absolute magnitudes in filters in array kf at z = 0
    vm=[]
    if bt.bolflux <=0:
        vm.append(-99.)
    else:
        vm.append(4.75-2.5*np.log10(bt.bolflux))       #    Bolometric magnitude
    for k in range(len(kf)):
        if kf[k] > 0:
            # Vega magnitude of model sed at redshift 0
            vm.append(VEGAmag(k,x,y,0.))
        else:
            # AB magnitude of model sed at redshift 0
            vm.append(ABmag(k,x,y,0.))
    # For the FUV, NUV and 1500 A bands, compute flux F(FUV), F(NUV), F(1500A) in Jansky
    aux = vm[len(vm)-1]
    vm[len(vm)-1] = flux_jansky(vm[0])
    vm.append(flux_jansky(vm[2]))
    vm.append(flux_jansky(aux))
    return vm

def ABmag(i,w,y,z):
    # Compute absolute AB system magnitude for SED (w,y) in filter i (via S. Charlot)
    # dl = 10 pc in units of Mpc
    dl = 1.e-5
    # This factor is SQRT(4*pi*(3.0856E24)^2/Lsun)
    AB0 = 5. * np.log10(1.7684e+08 * dl)
    # Compute flux F_mean
   #fmean = of_mean(i,w,y,z)
    fmean = qf_mean(i,w,y,z)
    # Compute AB absolute magnitude
    if fmean <= 0:
        ABmag = -99.
    else:
        ABmag = AB0 - 2.5*np.log10(fmean) - 48.6
    return ABmag

def VEGAmag(i,w,y,z):
    # Compute absolute VEGA system magnitude for SED (w,y) in filter i (via S. Charlot)
    fd = bt.ffd
    zp = fd[9]		# Zero points with respect to Vega
    # Zero point of ith filter
    f0 = zp[i]
    # Compute number of photons
   #fluxn = ofilter_n(i,w,y,z)
    fluxn = qfilter_n(i,w,y,z)
    # Compute Vega mag
    if fluxn <= 0:
        VEGAmag = -99.
    else:
        VEGAmag = f0 - 2.5*np.log10(fluxn)
    bt.filtern = fluxn
    return VEGAmag

def STmag(i,w,y,z):
    # Compute absolute STmag and corresponding flux for selected SED and filters
    # Compute flux F_meanl
    fmeanl = qf_meanl(i,w,y,z)
    # Compute ST mag
    if fmeanl <= 0:
        STmag = -99.
    else:
        STmag  = 16.232 - 21.10 - 2.5*np.log10(fmeanl)
    bt.fmeanl = fmeanl
    return STmag

def flux_jansky(ABmag):
    # Returns flux f_nu corresponding to a given ABmag in Jansky
    # 1 Jansky == 10^(−26) W/m²/Hz (MKS), or 10^(−23) ergs/s/cm²/Hz (CGS)
    # Ref:  https://pdn4kd.github.io/2020/10/28/janskyabmag.html
    if ABmag > -99:
        flux_jansky = 10.**(3.56 - 0.4*ABmag)
    else:
        flux_jansky = 0.
    return flux_jansky

def mass2light(mass,bmag,vmag,kmag):
    # Compute mass-to-light ratio in solar units. Improved definition (Feb. 27, 2004).
    # Sun absolute V, B and K magnitude
    vsun = 4.80
    bsun = 5.45
    ksun = 3.29
    # => (B-V)sun = 0.65, (V-K)sun = 1.51

    # M/L corresponding to 'mass'
    mlb = mass*10.**(0.4*(bmag-bsun))
    mlv = mass*10.**(0.4*(vmag-vsun))
    mlk = mass*10.**(0.4*(kmag-ksun))
    return mlb,mlv,mlk

def lx(x,y):
    # Compute X ray luminosity in ergs/sec between
    # 0.5 and 8 keV (i.e. between 1.5498 and 24.7969 Å)
    # Assumes X is in Angstroms and Y in Lo/Angstroms
    x1 = 1.5498 ; x2 = 24.7969 # ; lsun = 3.846E33
    if x[0] >= x2:
        lx = 0.
    else:
        xl = max(x1,x[0])
        lx = bs.trapz2(x,y,xl,x2)
    return lx

def ionphot(x,y):
    # Compute number of Lyman continuum photons in sed Y(X).
    # Compute number of Helium ionizing photons in sed Y(X).
    # Assumes X is in Angstroms and Y in ergs/sec/Angstroms (physical flux)

    # Physical constants
    c     = 2.997925E10
    h     = 6.6262E-27
    const = 1.0E-8/h/c
    sl    = 33. + np.log10(3.826)	# log of solar luminosity

    # Find Lyman limit in sed and number of photons
    wly   = 912.
    if x[0] >= wly:
        # clyman=0.
        # chelium=0.
        # che2=0.
        return 0., 0., 0.
    w=[]
    f=[]
    for i in range(len(x)):
        if x[i] < wly:
            w.append(x[i])
            f.append(x[i]*y[i])
        elif x[i] == wly:
            w.append(x[i])
            f.append(x[i]*y[i])
            break
        elif x[i] > wly:
            w.append(wly)
            fi=y[i-1] + (wly-x[i-1])*(y[i]-y[i-1])/(x[i]-x[i-1])
            f.append(wly*fi)
            break
    clyman = const*np.trapz(f,w)
    if clyman > 0:
        clyman = sl + np.log10(clyman)
    else:
        clyman = 0.

    # Find Helium limit in sed and number of photons
    # First electron, all photons shortward of whelium = 504.4 A
    wheI = 504.4
    if x[0] >= wheI:
        # chelium=0.
        # che2=0.
        return clyman, 0., 0.
    w=[]
    f=[]
    for i in range(len(x)):
        if x[i] < wheI:
            w.append(x[i])
            f.append(x[i]*y[i])
        elif x[i] == wheI:
            w.append(x[i])
            f.append(x[i]*y[i])
            break
        elif x[i] > wheI:
            w.append(wheI)
            fi=y[i-1] + (wheI-x[i-1])*(y[i]-y[i-1])/(x[i]-x[i-1])
            f.append(wheI*fi)
            break
    chelium = const*np.trapz(f,w)
    if chelium > 0:
        chelium = sl + np.log10(chelium)
    else:
        chelium = 0.

    # Second electron, all photons shortward of whelium = 228. A
    wheII = 228.
    if x[0] >= wheII:
        # che2=0.
        return clyman, chelium, 0.
    w=[]
    f=[]
    for i in range(len(x)):
        if x[i] < wheII:
            w.append(x[i])
            f.append(x[i]*y[i])
        elif x[i] == wheII:
            w.append(x[i])
            f.append(x[i]*y[i])
            break
        elif x[i] > wheII:
            w.append(wheII)
            fi=y[i-1] + (wheII-x[i-1])*(y[i]-y[i-1])/(x[i]-x[i-1])
            f.append(wheII*fi)
            break
    che2 = const*np.trapz(f,w)
    if che2 > 0:
        che2 = sl + np.log10(che2)
    else:
        che2 = 0.

    # Ratio de He II to H I ionizing photons
    cher = che2 - clyman
    ion=[]
    ion.append(clyman)
    ion.append(chelium)
    ion.append(che2)
    ion.append(cher)
    return ion

def filterid(k):
    # Returns ID for filter No. k

    # Filter log in file
    naml = os.environ.get('glxaux') + '/filterfrm.log'

    # Get filter id in log file
    k=int(k)
    with open(naml,'r') as f:
        if k == 0:
            # List ID for all filters
            l=-1
            print()
            print(' Log of filters in file: ',naml)
            print()
            for line in f:
                l+=1
                if l > 0:
                    print(str(l)+': '+line.rstrip())
        else:
            # Return ID for filter k
            l=0
            for line in f:
                l+=1
                if l==k+1:
                    break
    f.close()
    return line.rstrip()

def ABphot(w,y,z,dm,tl,tz,av,a,o,j):
    # Returns ABmag and FluxJansky for selected SED and filters
    # Define local variables
    u = []
    v = []
    ip = bt.ffd[0]
    for i in range(len(ip)):
        abmag = ABmag(i,w,y,z)
        jansk = flux_jansky(abmag)
        u.append(abmag)
        v.append(jansk)
    hd.printmag(a,z,dm,tl,tz,av,u,v,o,j)

def VEGAphot(w,y,z,dm,tl,tz,av,a,o,j):
    # Returns Vega magnitude for selected SED and filters
    # Define local variables
    u = []
    v = []
    ip = bt.ffd[0]
    for i in range(len(ip)):
        vgmag = VEGAmag(i,w,y,z)
        u.append(vgmag)
        v.append(bt.filtern)
    hd.printmag(a,z,dm,tl,tz,av,u,v,o,j)

def STphot(w,y,z,dm,tl,tz,av,a,o,j):
    # Returns STmag and corresponding flux for selected SED and filters
    # Define local variables
    u = []
    v = []
    ip = bt.ffd[0]
    for i in range(len(ip)):
        stmag = STmag(i,w,y,z)
        u.append(stmag)
        v.append(bt.fmeanl)
    hd.printmag(a,z,dm,tl,tz,av,u,v,o,j)

def vega_0p_n(g):
    # Returns zero point with respect to VEGA sed
    # Vega sed is first expressed in units of Lsun/A, as in the *.ised files
    # The A0 V stellar s.e.d. is read from file A0VSED the first time the function is used.
    # Uses number of photons inside filter (Fukugita et al. 1996)

    def A0VSED():
        # Return VEGA scaled sed
        # Vega apparent V magnitude
        vmag=0.03
        # Vega bolometric correction (old value)
        # bc=-0.1164
        # Vega bolometric correction (from Lang, Astrophysical Data, Table 9.3, p. 118)
        bc=-0.25
        # Vega bolometric magnitude
        vegabol=vmag+bc
        # Desired total flux corresponding to mbol = vegabol in units of Lsun
        stot=10.**(-0.4*(vegabol-4.75))
        # Read A0VSED file
        a0 = bf.nfile(os.environ.get('A0VSED'),1)
        x0 = a0['col1']
        y0 = a0['col2']
        # Compute total flux below Vega sed
        # tot=trapz1(xa0v,ya0v,jread)
        tot = np.trapz(y0,x0)
        # Scale sed to desired stot, final sed in units of Lsun/A as *.ised files
        scl=stot/tot
        y0 = scl*y0
        return x0,y0

    # Read A0VSED file
    x0,y0 = A0VSED()

    # Read and organize filters in array g
    bt.ffd = qfilters(x0,g,0.,0.,0)

    # Compute zeropoints for all filters in ip
    # vega_0p_n=2.5*alog10(filter_n(mf,xa0v,yaux5,jread,0.,kerr))
    u = []
    ip = bt.ffd[0]
    for i in range(len(ip)):
        filtern = ofilter_n(i,x0,y0,0.)
        vgzp = 2.5*np.log10(filtern)
        u.append(vgzp)
    return u

def frm2fits():

    # Build fits file with filter response functions

    # Read filter file in text format and creates a fits file
    name = os.environ.get('glxaux') + '/filterfrm.res'
    namo = os.environ.get('glxaux') + '/filterfrm.out'
    namh = os.environ.get('glxaux') + '/filterfrm.hdr'
    namz = os.environ.get('glxaux') + '/filterfrm.vzp'
    naml = os.environ.get('glxaux') + '/filterfrm.log'
    namf = os.environ.get('glxaux') + '/FILTERRES.fits'
    print('Reading file: ' + name)
    print('Writing file: ' + namf)

    # Init parameters
    bt.iread = True

    # Write response functions to text file with no header information
    c = 0
    l = 0
    k = 0
    n = []
    h = []
    o = open(namo,'w')
    o.write('# Wavelength  Response\n')
    with open(name,'r') as f:
        # Read each file in the input list
        for line in f:
            l+=1
            if '#' in line:
                c+=1
                q = True
                line = line.replace('# ','')
                h.append(line)
            else:
                k+=1
                o.write(line)
                if q:
                    n.append(k)
                    q = False
    f.close()
    o.close()
    # print (c,'filters in file')
    # print (l,'lines in file')
    # print (k,'lines in out file')
    # print (len(n))
    # print(n)

    # Write filter headers and filter log to text files
    t = open(namh,'w')
    d = open(naml,'w')
    t.write('#  N         R1          R2    Zeropoint\n')
    d.write('#  Filter_ID\n')
    for i in range(len(n)):
        n1 = n[i]
        if i < len(n)-1:
            n2 = n[i+1]-1
        else:
            n2 = k
        # Zero point = 0. Later on will be computed.
        fz = 0.
        l = f'{i+1:4d} {n1:10d}  {n2:10d} {fz:12.6f}' + '\n'
        t.write(l)
        d.write(h[i])
    t.close()
    d.close()

    # Create fits file with header information and filter transmissions
    t = Table.read(namh, format='ascii')
    t.write(namf, overwrite=True)
    t = Table.read(namo, format='ascii')
    t.write(namf, append=True)

    # Compute zeropoints in the Vega magnitude system for the n filters using
    # the fits file created above. Add zeropoints to new header file
    o = open(namz,'w')
    o.write('#  N         R1          R2    Zeropoint\n')
    h = bf.nfile(namh,0)
    g=[]
    for i in range(c):
        g.append(i+1)
    f0=vega_0p_n(g)
    for i in range(c):
        j1 = h[i][0]
        j2 = h[i][1]
        j3 = h[i][2]
        fz = f0[i]
        l = f'{j1:4d} {j2:10d}  {j3:10d} {fz:12.6f}' + '\n'
        o.write(l)
    o.close()

    # Create new fits file
    t = Table.read(namz, format='ascii')
    t.write(namf, overwrite=True)
    t = Table.read(namo, format='ascii')
    t.write(namf, append=True)

    # List filters in file:
    filterid(0)

    # Delete temporary files
    os.system('\\rm -f ' + namo + ' ' + namh + ' ' + namz)


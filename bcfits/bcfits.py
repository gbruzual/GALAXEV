import os
import sys
import warnings
import numpy as np
from   scipy.io import FortranFile
from   astropy.table import Table
from   astropy.io import fits
from   astropy.utils.exceptions import AstropyWarning
import astropy.io.ascii as ascii
import builtins        as bt
import basics.basics   as bs
import bchead.bchead   as hd

def nfile(file,i):
    # Read ascii file
    if i==0:
        print ('Reading file: ',file)
    d = ascii.read(file)
    return d

def bcfits(file):
    # Reads fits table with BC models if different from last read
    if file != bt.lfits:
        bt.bc = read_bcfits(file)
    return bt.bc

def read_bcfits(file):
    # Reads fits table with BC models

    # Read file in fits table format
    file = bs.fcheck(file)
    with fits.open(file) as hdul:
        warnings.simplefilter('ignore', category=AstropyWarning)
        # Use Table to read the fits table
        f = Table.read(hdul,hdu=1)	# hdu = 1 => SED's (luminosity vs. wavelength)
        p = Table.read(hdul,hdu=2)	# hdu = 2 => galaxy physical properties
        m = Table.read(hdul,hdu=3)	# hdu = 3 => photometric magnitude in different bands
        d = Table.read(hdul,hdu=4)	# hdu = 4 => line spectral indices
        t = Table.read(hdul,hdu=5)	# hdu = 5 => time scale for spectral evolution (221 steps)
        s = Table.read(hdul,hdu=6)	# hdu = 6 => miscellaneous data concerning IMF and SFR
        v = t                           # hdu = 5 => time scale, includes flux weighted ages
        t = t['age-yr']			# time scale without pre MS evolution (agrees with t for BC03 and C&B models)
        w = f['Wavelength']		# wavelength array
        a = 10.**m['logage']		# time scale without pre MS evolution (agrees with t for BC03 and C&B models)
        if a[0] < 10.:
            a[0]=0.			# change first time step to a = 0. It is written as 1 yr in the color files.
        e = hdul[0].header		# ascii table header
        h = hdul[1].header		# SED column header
        for i in range(len(t)+1):
            h[i] = h['TTYPE' + str(i+1)]
   #return w,f,t,e,s,m,p,d,v,a,h
    bt.bc = w,f,t,e,s,m,p,d,v,a,h
    bt.lfits = bs.lnam(file)
    return bt.bc

def age_yr(age):
    # Checks if age was entered in yr or Gyr. Returns age in yr.
    if age <= 20.:
        # age was entered in Gyr
        age = age*1.E9
    return age

def ft(age,f,t):
    # Return flux in record corresponding to age t in yr
    age = age_yr(age)
    i = np.abs(np.array(t) - age).argmin()    # index of element of array t closest in value to age a
    y = f[i+1][:]
    return y

def fi(age,f,t):
    # Return flux in record corresponding to age t interpolating if necessary
    age = age_yr(age)
    i1 = np.searchsorted(t,age)
    if t[i1] == age:
        # Record of desired age exists
        y=ft(age,f,t)
        # print('->',age,i1,t[i1])
    else:
        # Interpolate record
        i2=i1-1
        if age > 0 and t[i1] > 0 and t[i2] > 0:
            a1=np.log10(t[i2]/age)/np.log10(t[i2]/t[i1])
        else:
            a1=(t[i2]-age)/(t[i2]-t[i1])
        a2=1.-a1
        y1=ft(t[i1],f,t)
        y2=ft(t[i2],f,t)
        y =a1*y1+a2*y2
        # print('-->',age,i2,t[i2],i1,t[i1],a2,a1)
    return y

def pt(age,p,m,t,i):
    # Return value of i-th physical property at age = t
    age = age_yr(age)
    j = np.abs(np.array(t) - age).argmin()    # index of element of array t closest in value to age
    if i==1:
        y = 10.**(-0.4*(m[j][1]-4.75))			# bolometric flux
    else:
        y = p[j][i]
    return y

def pi(age,p,m,t,i):
    # Return value of i-th physical property at age = t, interpolating if necessary
    age = age_yr(age)
    i1 = np.searchsorted(t,age)
    if t[i1] == age:
        # Record of desired age exists
        if i==1:
            d = np.array(m[i1][1])
            y = float(d)
            y = 10.**(-0.4*(y-4.75))			# bolometric flux
        else:
            y = p[i1][i]
    else:
        # Interpolate record
        i2=i1-1
        if age > 0 and t[i1] > 0 and t[i2] > 0:
            a1=np.log10(t[i2]/age)/np.log10(t[i2]/t[i1])
        else:
            a1=(t[i2]-age)/(t[i2]-t[i1])
        a2=1.-a1
        y1=pt(t[i1],p,m,t,i)
        y2=pt(t[i2],p,m,t,i)
        y =a1*y1+a2*y2
    return y

def bcnew(fo):
    # Interpolate BC/CB model to the required metallicity
    global b1, a1, b2, a2, bp1, bp2
    global x, t, y1,t1,m1,p1,s1, y2,t2,m2,p2,s2

    # Check if interpolation in metallicity has been requested
    if ',' in fo:
        i  = fo.find(',')
        zs = fo[i+1:]
        zx = float(zs)
        fo = fo[:i]

        # Available metallicity values
        z03 = [0.0001, 0.0004,  0.004, 0.008,  0.02,  0.05,  0.10]   						                    # BC03/BC07
        a03 = ['z0001','z0004', 'z004','z008', 'z02', 'z05', 'z10']   						                    # BC03/BC07
        z19 = [0.0000, 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.004, 0.006, 0.008, 0.010, 0.014, 0.017, 0.020, 0.030, 0.040, 0.060]  # CB19
        a19 = ['z0000','z0001','z0002','z0005','z001','z002','z004','z006','z008','z010','z014','z017','z020','z030','z040','z060'] # CB19

        # MockGal052      1.523     0.233     0.806     1.609     0.177     0.465     0.286     1.329     0.743     0.324      393.

        # Find requested model
        if 'bc2003' in fo or 'cb2007' in fo:
            b03  = True
            b19  = False
            zmod = z03
            amod = a03
            zsun = 0.02
        elif 'bc2019' in fo or 'cb2019' in fo:
            b03  = False
            b19  = True
            zmod = z19
            amod = a19
            zsun = 0.017
        else:
            print(' Unknown type of model:',fo)

        # Express zx in units of Zsun
        if zx < 0:
            # If zx < 0, zx is in units of Zsun
            zx = -zx*zsun
        zs = f'{zx:6.4f}'
        zs = zs.replace('0.','z')

        # Find metallicities bracketing zx and build corresponding file names
        i  = fo.find('_')
        j  = fo.find('_',i+1)
        ix = np.searchsorted(zmod,zx)
        if ix == len(zmod):
            fo = fo[:i+1] + amod[ix-1] + fo[j:]
        elif ix == 0:
            fo = fo[:i+1] + amod[ix]   + fo[j:]
        else:
            f1 = fo[:i+1] + amod[ix-1] + fo[j:]
            f2 = fo[:i+1] + amod[ix]   + fo[j:]
            fo = fo[:i+1] + zs         + fo[j:]
            f1 = bs.fcheck(f1)
            f2 = bs.fcheck(f2)
           #x,y1,t1,e1,s1,m1,*nouse = read_bcfits(f1)
           #x,y2,t2,e2,s2,m2,*nouse = read_bcfits(f2)
            x,y1,t1,e1,s1,m1,*nouse = bcfits(f1)
            x,y2,t2,e2,s2,m2,*nouse = bcfits(f2)
            t  = t1
            # Use addbursts to perform interpolation in log(z)
            z1 = np.log10(zmod[ix-1])
            z2 = np.log10(zmod[ix])
            zo = np.log10(zx)
            a1 = (z2-zo)/(z2-z1)
            a2 = 1. - a1
            # Define burst parameters
            b1,b2 = 0,0
            bp1 = [b1,a1,f1]
            bp2 = [b2,a2,f2]

            # Add 2 bursts
            addbursts('NEW Z MODEL',fo)

    # Read BC/CB model in fits file
    fo = bs.fcheck(fo)
   #x,y,t,e,s,m,p,*nouse = read_bcfits(fo)
    x,y,t,e,s,m,p,*nouse = bcfits(fo)
    hd.multhead(0,fo,t,x,e)
    return x,y,t,e,s,m,p

def read_bcised(ifile):

    # Read *.ised binary file written by the galaxev fortran code

    # Comments by GBA 02.02.2020 (Madrid):
    #     The read_record function in FortranFile does not accept combining integer and float values in a single read
    #     A practical solution is to open the file once to read nt, nw, and ni as integer numbers, and then
    #     open the file again to read:
    #         t = time scale (nt steps)
    #         w = wavelength (nw points)
    #         f = sed fluxes (nw points, nt times)
    #         h = line index (ni points, nt times)
    #         p = phys prop  (nt points, 16 properties)

    # Open file to read 3 first records as integer
    with FortranFile(ifile) as file:
        # Read 1st record and get number nt of time steps
        d = file.read_record(dtype=np.int32)
        nt = d[0]
        # Read 2nd record and get number nw of wavelength points
        d = file.read_record(dtype=np.int32)
        nw = d[0]
        # Read 3rd record and get number ni of line index fluxes written after sed
        d = file.read_record(dtype=np.int32)
        ni = d[nw+1]
        # print(nt,' time steps')
        # print(nw,' wavelength points')
        # print(ni,' line index fluxes')

    # Open file again to read all record as float
    with FortranFile(ifile) as file:
        # Read 1st record with time scale
        d = file.read_record(dtype=np.float32)
        t = d[1:nt+1]
        h = []
        h.append(t[1:nt+1])			# flux
        p=[]
        p.append(t[1:nt+1])			# flux
        # Read 2nd record with wavelength scale
        d = file.read_record(dtype=np.float32)
        w = d[1:nw+1]
        f = []
        f.append(w[1:nw+1])			# flux
        # Read nt records with flux and line index evolution
        for i in range(nt):
            d = file.read_record(dtype=np.float32)
            f.append(d[1:nw+1])			# flux
            if nw != 56:
                h.append(d[nw+2:])		# line indices (not for ised files in JPAS filter system)
        # Read 16 records with physical properties for this model
        if nw != 56:
            for i in range(16):
                d = file.read_record(dtype=np.float32)
                p.append(d[1:])			# physical properties
        # Build array with number of time step
        a = []
        for i in range(len(t)+1):
            a.append(i)
    # Report file content
    hd.pid(ifile,w,t)
   #print()
   #print("In file " + bs.bp(bs.lnam(ifile)) + " there are " + str('{:d}'.format(len(t))) + " galaxy SED's, ranging from " + str('{:.1E}'.format(t[0])) +
   #      " to " + str('{:.1E}'.format(t[len(t)-1])) + " yr")
   #print("Each SED covers the wavelength range from " + str('{:.1f}'.format(w[0])) + " to " + str('{:.1E}'.format(w[len(w)-1])) +
   #      " A in " + str('{:d}'.format(len(w))) + " steps")

    # Return arrays
    # return w,f,t,h,p		# uncomment if h and p are needed
    return w,f,t,a

def read_bcased(ifile):

    # Read *.ased ascii file written by the galaxevpl fortran code

    # Read time scale from header
    k = 0
    print()
    print('Reading text file: ' + bs.bp(bs.lnam(ifile)) + ' will take a few seconds')
    with open(ifile, 'r') as f:
        for line in f:
            k = k+1
            h = line
            if k == 6:
                h = h.replace('# Lambda(A)',' 0.00E+00  ')
                h = h.split()
                h = np.array(h,dtype=np.float32)
                t = [0]
                t = h[1:]
                break

    # Read seds at all ages
    f = ascii.read(ifile)
    w = f[0][:]

    # Build array with number of time step
    a = []
    for i in range(len(t)+1):
        a.append(i)

    # Report file content
    hd.pid(ifile,w,t)
   #print("In file " + bs.bp(bs.lnam(ifile)) + " there are " + str('{:d}'.format(len(t))) + " galaxy SED's, ranging from " + str('{:.1E}'.format(t[0])) +
   #      " to " + str('{:.1E}'.format(t[len(t)-1])) + " yr")
   #print("Each SED covers the wavelength range from " + str('{:.1f}'.format(w[0])) + " to " + str('{:.1E}'.format(w[len(w)-1])) +
   #      " A in " + str('{:d}'.format(len(w))) + " steps")
    return w,f,t,a

def wrfits(f,k):
    # Write fits file with model results
    global io

    # Add IMF miscelaneous data to fits file
    for i in range(len(s1)):
        sx = bt.s1[:][i]
        sw =[]
        for j in range(10):
            sw.append(sx[j])
        bt.td6.add_row(sw)
    iop = sw[2]
    if k==1:
        # For add burst code store IMF for second model
        for i in range(len(s2)):
            sx = bt.s2[:][i]
            sw =[]
            for j in range(10):
                sw.append(sx[j])
            bt.td6.add_row(sx)
        # Store burst parameters
        bt.td7.add_row(bt.bp1)
        bt.td7.add_row(bt.bp2)

    # Write fits file
    print()
    print()
    f = f.replace('.fits','')
    g = f
    o = os.environ.get('glxtmp') + '/'
    o = o.replace('//','/')
    f = o + f + '.fits'
    t = open(o + '__lastfits__.txt','w')
    t.write(f)
    t.close()
    o = os.environ.get('glxout') + '/'
    o = o.replace('//','/')
    f = o + g + '.fits'
    print(' Creating file(s): ' + f)
    bt.td1.write(f, overwrite=True)
    bt.td2.write(f, append=True)
    bt.td3.write(f, append=True)
    bt.td4.write(f, append=True)
    bt.td5.write(f, append=True)
    bt.td6.write(f, append=True)
    if k==1:
        bt.io = 5
        bt.td7.write(f, append=True)
    # Update SED headers in fits file
    rbldfitshdr(f,iop)
    print()
    if bt.io == 9:
        h = o + 'py.500'
        g = o + g + '.sfr'
        os.system('\\mv -f ' + h + '  ' + g)
        print(' SFR saved in table:',g)
        print()

def rbldfitshdr(name,iop):
    # Rebuilds time scale header for SEDs in BC/CB model in fits file
    import os

    # Read fits table with BC models. Recover time steps
    hdul   = fits.open(name)
    t = Table.read(hdul,hdu=5)	# hdu = 5 => time scale
    t = t['age-yr']

    # Build fits compatible header for BINTABLE1 (model age)	
    sechdr = hdul[1].header
    for i in range(len(t)+1):
        if i==0:
            s = 'Wavelength'
            u = 'A  (Wavelength)'
        else:
            c = t[i-1]
            if c < 10.E9:
                s = str('{:.7E}'.format(c))
            else:
                s = str('{:.6f}'.format(c*1.E-9)) + "E9"
            s = s.replace("E+1", "E")
            u = 'Lo/A  (SED at t = ' + s + ' yr)'
            s = s.replace(".", "p")
            s = 't' + s
        sechdr['TTYPE' + str(i+1)]  = s 	# label for column i+1
        sechdr['TUNIT' + str(i+1)]  = u		# units for column i+1

    # Update headers
    prihdr = hdul[0].header
    prihdr['stilts0']  = ''
    prihdr['stilts']   = 'Fits file created by STILTS v3.1-1, including header up to NTABLE'
    prihdr['modlhd0']  = ''
    prihdr['modlhdr']  = '---------------- CB2016 model header follows below ----------------'
    prihdr['biblio00']  = ''
    prihdr['biblio01']  = 'BIBLIOGRAPHY:'
    prihdr['biblio02']  = ''
    prihdr['biblio03']  = 'EVOLUTIONARY TRACKS:'
    prihdr['biblio04']  = '             PARSEC: Bressan, A., et al. 2012, MNRAS, 427'
    prihdr['biblio05']  = '                     Chen, Y., et al. 2015, MNRAS, 452, 1068'
    prihdr['biblio06']  = '             TP-AGB: Marigo, P. et al. 2013, MNRAS, 434, 488'
    prihdr['biblio07']  = ''
    prihdr['biblio08']  = 'STELLAR SPECTRA:'
    prihdr['biblio09']  = '      Visible range: MILES stellar library'
    prihdr['biblio10']  = '                     Sanchez-Blazquez et al. 2006, MNRAS, 371, 703'
    prihdr['biblio11']  = '                     Falcon-Barroso, J., et al. 2011, A&A, 532, 95'
    prihdr['biblio12']  = '                     Prugniel, Ph. 2011 A&A, 531, 165'
    prihdr['biblio13']  = '  Supplemented with:'
    prihdr['biblio14']  = '             Stelib: redward of Miles reddest point (7350 A)'
    prihdr['biblio15']  = '                     Le Borgne, J.-F. et al. 2003, A&A, 402, 433'
    prihdr['biblio16']  = '          BaSeL 3.1: redward of Stelib reddest point (8750 A)'
    prihdr['biblio17']  = '                     Westera P. et al. 2002, A&A, 381, 524'
    prihdr['biblio18']  = '            O-stars: Tlusty models'
    prihdr['biblio19']  = '                     Lanz, T. & Hubeny, I. 2003, ApJS, 146, 417'
    prihdr['biblio20']  = '            B-stars: Tlusty models'
    prihdr['biblio21']  = '                     Lanz, T. & Hubeny, I. 2007, ApJS, 169, 83'
    prihdr['biblio22']  = '   Massive MS stars: WM-Basic models'
    prihdr['biblio23']  = '                     Leitherer, C., et al. 2010, ApJS, 189, 309'
    prihdr['biblio24']  = '           WR stars: PoWR models'
    prihdr['biblio25']  = '                     Hamann, W-R & Gafener, G. 2004, A&A, 427, 697'
    prihdr['biblio26']  = '            A-stars: Martins et al. models'
    prihdr['biblio27']  = '                     Martins, L.P. et al. 2005, MNRAS, 358, 49'
    prihdr['biblio28']  = '   Lower mass stars: UV-Blue models'
    prihdr['biblio29']  = '                     Rodriguez-Merino, LH, et al. 2005,ApJ,626,411'
    prihdr['biblio30']  = '               CSPN: Rauch models'
    prihdr['biblio31']  = '                     Rauch, T. 2002, RevMexAstronAstrof CS, 12, 150'
    prihdr['biblio32']  = '       TP-AGB stars: Aringer models for C stars'
    prihdr['biblio33']  = '                     Aringer et al. 2009, A&A, 503, 913'
    prihdr['biblio34']  = '                     IRTF library'
    prihdr['biblio35']  = '                     Rayner et al. 2009, ApJS, 185, 289'
    prihdr['biblio36']  = '                     Dusty code'
    prihdr['biblio37']  = '                     Ivezic Z. & Elitzur M., 1997, MNRAS, 287, 799'
    prihdr['biblio38']  = '                     Gonzalez-L. R. et al. 2010, MNRAS, 403, 1213'
    prihdr['biblio39']  = ''
    prihdr['biblio40']  = 'LINE INDICES:'
    prihdr['biblio41']  = '       Lick indices: Worthey, G., et al. 1994, ApJS, 94, 687'
    prihdr['biblio42']  = '                     Trager, S.C., et al. 1998, ApJS, 116, 1'
    prihdr['biblio43']  = '      Other indices: Worthey G & Ottaviani DL 1997, ApJS, 111, 377'
    prihdr['biblio44']  = '                     Brodie, J., & Hanes, D.A. 1986, ApJ, 300, 258'
    prihdr['biblio45']  = '                     Diaz, A. & Terlevich**2 1989, MNRAS, 239, 325'
    prihdr['biblio46']  = '                     Gorgas J, Cardiel N et al. 1999, A&AS, 139, 29'
    prihdr['biblio47']  = '                     Marcillac, D., et al. 2006, A&A, 458, 369'
    prihdr['biblio48']  = '         UV indices: Fanelli, M. N., et al. 1987, ApJ, 321, 768'
    prihdr['biblio49']  = '                     Fanelli, M. N., et al. 1990, ApJ, 364, 272'
    prihdr['biblio50']  = '                     Fanelli, M. N., et al. 1992, ApJS, 82, 197'
    prihdr['biblio51']  = '                     Chavez, M. et al. 2009, ApJ, 700, 694'
    prihdr['biblio52']  = '                     Maraston, C. et al. 2009, A&A, 493, 425'
    prihdr['biblio53']  = ''
    prihdr['biblio54']  = 'Fuel Consumption Theorem: Renzini & Buzzoni 1986, in'
    prihdr['biblio55']  = '  Spectral Evolution of Galaxies, eds. C. Chiosi & A. Renzini, p195'
    prihdr['biblio56']  = ''
    prihdr['biblio57']  = '-------------------------------------------------------------------'
    prihdr['modlhd1']  =   ''
    prihdr['modlid']   = 'CB2016 - Simple Stellar Population Synthesis Model'
    prihdr['modlref']  = 'Charlot & Bruzual (2016, MNRAS, in preparation)'
    prihdr['modlfile'] = str(name)

    # BC2003 models, Padova 1994 tracks + VWPAGB
    if str(iop) == '-22':
        prihdr['modltrks'] = ('0.7696, 0.2303, 0.0001   ','X, Y, Z, Padova 1994 + TP-AGB + VWPAGB')
    elif str(iop) == '-32':
        prihdr['modltrks'] = ('0.7686, 0.2310, 0.0004   ','X, Y, Z, Padova 1994 + TP-AGB + VWPAGB')
    elif str(iop) == '-42':
        prihdr['modltrks'] = ('0.7560, 0.2400, 0.0040   ','X, Y, Z, Padova 1994 + TP-AGB + VWPAGB')
    elif str(iop) == '-52':
        prihdr['modltrks'] = ('0.7420, 0.2500, 0.0080   ','X, Y, Z, Padova 1994 + TP-AGB + VWPAGB')
    elif str(iop) == '-62':
        prihdr['modltrks'] = ('0.7000, 0.2800, 0.0200   ','X, Y, Z, Padova 1994 + TP-AGB + VWPAGB')
    elif str(iop) == '-72':
        prihdr['modltrks'] = ('0.5980, 0.3520, 0.0500   ','X, Y, Z, Padova 1994 + TP-AGB + VWPAGB')
    elif str(iop) == '-82':
        prihdr['modltrks'] = ('0.4250, 0.4750, 0.1000   ','X, Y, Z, Padova 1994 + TP-AGB + VWPAGB')

    # BC2003 models, Padova 2000 tracks + VWPAGB
    elif str(iop) == '-122':
        prihdr['modltrks'] = ('0.7696, 0.2300, 0.0004   ','X, Y, Z, Padova 2000 + TP-AGB + VWPAGB')
    elif str(iop) == '-132':
        prihdr['modltrks'] = ('0.7690, 0.2300, 0.0010   ','X, Y, Z, Padova 2000 + TP-AGB + VWPAGB')
    elif str(iop) == '-142':
        prihdr['modltrks'] = ('0.7560, 0.2400, 0.0040   ','X, Y, Z, Padova 2000 + TP-AGB + VWPAGB')
    elif str(iop) == '-152':
        prihdr['modltrks'] = ('0.7420, 0.2500, 0.0080   ','X, Y, Z, Padova 2000 + TP-AGB + VWPAGB')
    elif str(iop) == '-162':
        prihdr['modltrks'] = ('0.7080, 0.2730, 0.0190   ','X, Y, Z, Padova 2000 + TP-AGB + VWPAGB')
    elif str(iop) == '-172':
        prihdr['modltrks'] = ('0.6700, 0.3000, 0.0300   ','X, Y, Z, Padova 2000 + TP-AGB + VWPAGB')

    # CB2007 models, Padova 1994 tracks + VWPAGB + enhanced TP-AGB
    elif str(iop) == '-92':
        prihdr['modltrks'] = ('0.7696, 0.2303, 0.0001   ','X, Y, Z, Padova 1994 + TP-AGB + VWPAGB')
    elif str(iop) == '-93':
        prihdr['modltrks'] = ('0.7686, 0.2310, 0.0004   ','X, Y, Z, Padova 1994 + TP-AGB + VWPAGB')
    elif str(iop) == '-94':
        prihdr['modltrks'] = ('0.7560, 0.2400, 0.0040   ','X, Y, Z, Padova 1994 + TP-AGB + VWPAGB')
    elif str(iop) == '-95':
        prihdr['modltrks'] = ('0.7420, 0.2500, 0.0080   ','X, Y, Z, Padova 1994 + TP-AGB + VWPAGB')
    elif str(iop) == '-96':
        prihdr['modltrks'] = ('0.7000, 0.2800, 0.0200   ','X, Y, Z, Padova 1994 + TP-AGB + VWPAGB')
    elif str(iop) == '-97':
        prihdr['modltrks'] = ('0.5980, 0.3520, 0.0500   ','X, Y, Z, Padova 1994 + TP-AGB + VWPAGB')
    elif str(iop) == '-98':
        prihdr['modltrks'] = ('0.4250, 0.4750, 0.1000   ','X, Y, Z, Padova 1994 + TP-AGB + VWPAGB')

    # BC2019 models, PARSEC tracks + VWPAGB
    elif str(iop) == '-400':
        prihdr['modltrks'] = ('0.7700, 0.2300, 0.0000   ','X, Y, Z, PARSEC + TP-AGB + VWPAGB')
    elif str(iop) == '-401':
        prihdr['modltrks'] = ('0.7509, 0.2490, 0.0001   ','X, Y, Z, PARSEC + TP-AGB + VWPAGB')
    elif str(iop) == '-402':
        prihdr['modltrks'] = ('0.7508, 0.2490, 0.0002   ','X, Y, Z, PARSEC + TP-AGB + VWPAGB')
    elif str(iop) == '-403':
        prihdr['modltrks'] = ('0.7505, 0.2490, 0.0005   ','X, Y, Z, PARSEC + TP-AGB + VWPAGB')
    elif str(iop) == '-404':
        prihdr['modltrks'] = ('0.7490, 0.2500, 0.0010   ','X, Y, Z, PARSEC + TP-AGB + VWPAGB')
    elif str(iop) == '-405':
        prihdr['modltrks'] = ('0.7460, 0.2520, 0.0020   ','X, Y, Z, PARSEC + TP-AGB + VWPAGB')
    elif str(iop) == '-406':
        prihdr['modltrks'] = ('0.7400, 0.2560, 0.0040   ','X, Y, Z, PARSEC + TP-AGB + VWPAGB')
    elif str(iop) == '-407':
        prihdr['modltrks'] = ('0.7350, 0.2590, 0.0060   ','X, Y, Z, PARSEC + TP-AGB + VWPAGB')
    elif str(iop) == '-408':
        prihdr['modltrks'] = ('0.7290, 0.2630, 0.0080   ','X, Y, Z, PARSEC + TP-AGB + VWPAGB')
    elif str(iop) == '-409':
        prihdr['modltrks'] = ('0.7230, 0.2670, 0.0100   ','X, Y, Z, PARSEC + TP-AGB + VWPAGB')
    elif str(iop) == '-410':
        prihdr['modltrks'] = ('0.7130, 0.2730, 0.0140   ','X, Y, Z, PARSEC + TP-AGB + VWPAGB')
    elif str(iop) == '-411':
        prihdr['modltrks'] = ('0.7040, 0.2790, 0.0170   ','X, Y, Z, PARSEC + TP-AGB + VWPAGB')
    elif str(iop) == '-412':
        prihdr['modltrks'] = ('0.6960, 0.2840, 0.0200   ','X, Y, Z, PARSEC + TP-AGB + VWPAGB')
    elif str(iop) == '-413':
        prihdr['modltrks'] = ('0.6680, 0.3020, 0.0300   ','X, Y, Z, PARSEC + TP-AGB + VWPAGB')
    elif str(iop) == '-414':
        prihdr['modltrks'] = ('0.6390, 0.3210, 0.0400   ','X, Y, Z, PARSEC + TP-AGB + VWPAGB')
    elif str(iop) == '-416':
        prihdr['modltrks'] = ('0.5840, 0.3560, 0.0600   ','X, Y, Z, PARSEC + TP-AGB + VWPAGB')

    # CB2019 models, PARSEC tracks + MBPAGB
    elif str(iop) == '-500':
        prihdr['modltrks'] = ('0.7700, 0.2300, 0.0000   ','X, Y, Z, PARSEC + TP-AGB + MBPAGB')
    elif str(iop) == '-501':
        prihdr['modltrks'] = ('0.7509, 0.2490, 0.0001   ','X, Y, Z, PARSEC + TP-AGB + MBPAGB')
    elif str(iop) == '-502':
        prihdr['modltrks'] = ('0.7508, 0.2490, 0.0002   ','X, Y, Z, PARSEC + TP-AGB + MBPAGB')
    elif str(iop) == '-503':
        prihdr['modltrks'] = ('0.7505, 0.2490, 0.0005   ','X, Y, Z, PARSEC + TP-AGB + MBPAGB')
    elif str(iop) == '-504':
        prihdr['modltrks'] = ('0.7490, 0.2500, 0.0010   ','X, Y, Z, PARSEC + TP-AGB + MBPAGB')
    elif str(iop) == '-505':
        prihdr['modltrks'] = ('0.7460, 0.2520, 0.0020   ','X, Y, Z, PARSEC + TP-AGB + MBPAGB')
    elif str(iop) == '-506':
        prihdr['modltrks'] = ('0.7400, 0.2560, 0.0040   ','X, Y, Z, PARSEC + TP-AGB + MBPAGB')
    elif str(iop) == '-507':
        prihdr['modltrks'] = ('0.7350, 0.2590, 0.0060   ','X, Y, Z, PARSEC + TP-AGB + MBPAGB')
    elif str(iop) == '-508':
        prihdr['modltrks'] = ('0.7290, 0.2630, 0.0080   ','X, Y, Z, PARSEC + TP-AGB + MBPAGB')
    elif str(iop) == '-509':
        prihdr['modltrks'] = ('0.7230, 0.2670, 0.0100   ','X, Y, Z, PARSEC + TP-AGB + MBPAGB')
    elif str(iop) == '-510':
        prihdr['modltrks'] = ('0.7130, 0.2730, 0.0140   ','X, Y, Z, PARSEC + TP-AGB + MBPAGB')
    elif str(iop) == '-511':
        prihdr['modltrks'] = ('0.7040, 0.2790, 0.0170   ','X, Y, Z, PARSEC + TP-AGB + MBPAGB')
    elif str(iop) == '-512':
        prihdr['modltrks'] = ('0.6960, 0.2840, 0.0200   ','X, Y, Z, PARSEC + TP-AGB + MBPAGB')
    elif str(iop) == '-513':
        prihdr['modltrks'] = ('0.6680, 0.3020, 0.0300   ','X, Y, Z, PARSEC + TP-AGB + MBPAGB')
    elif str(iop) == '-514':
        prihdr['modltrks'] = ('0.6390, 0.3210, 0.0400   ','X, Y, Z, PARSEC + TP-AGB + MBPAGB')
    elif str(iop) == '-516':
        prihdr['modltrks'] = ('0.5840, 0.3560, 0.0600   ','X, Y, Z, PARSEC + TP-AGB + MBPAGB')

    # BC2021 models, SISSA 2020 tracks + VWPAGB
    elif str(iop) == '-600':
        prihdr['modltrks'] = ('0.7700, 0.2300, 0.0000   ','X, Y, Z, SISSA 2020 + TP-AGB + VWPAGB')
    elif str(iop) == '-601':
        prihdr['modltrks'] = ('0.7509, 0.2490, 0.0001   ','X, Y, Z, SISSA 2020 + TP-AGB + VWPAGB')
    elif str(iop) == '-602':
        prihdr['modltrks'] = ('0.7508, 0.2490, 0.0002   ','X, Y, Z, SISSA 2020 + TP-AGB + VWPAGB')
    elif str(iop) == '-603':
        prihdr['modltrks'] = ('0.7505, 0.2490, 0.0005   ','X, Y, Z, SISSA 2020 + TP-AGB + VWPAGB')
    elif str(iop) == '-604':
        prihdr['modltrks'] = ('0.7490, 0.2500, 0.0010   ','X, Y, Z, SISSA 2020 + TP-AGB + VWPAGB')
    elif str(iop) == '-605':
        prihdr['modltrks'] = ('0.7460, 0.2520, 0.0020   ','X, Y, Z, SISSA 2020 + TP-AGB + VWPAGB')
    elif str(iop) == '-606':
        prihdr['modltrks'] = ('0.7400, 0.2560, 0.0040   ','X, Y, Z, SISSA 2020 + TP-AGB + VWPAGB')
    elif str(iop) == '-607':
        prihdr['modltrks'] = ('0.7350, 0.2590, 0.0060   ','X, Y, Z, SISSA 2020 + TP-AGB + VWPAGB')
    elif str(iop) == '-608':
        prihdr['modltrks'] = ('0.7290, 0.2630, 0.0080   ','X, Y, Z, SISSA 2020 + TP-AGB + VWPAGB')
    elif str(iop) == '-609':
        prihdr['modltrks'] = ('0.7230, 0.2670, 0.0100   ','X, Y, Z, SISSA 2020 + TP-AGB + VWPAGB')
    elif str(iop) == '-610':
        prihdr['modltrks'] = ('0.7130, 0.2730, 0.0140   ','X, Y, Z, SISSA 2020 + TP-AGB + VWPAGB')
    elif str(iop) == '-611':
        prihdr['modltrks'] = ('0.7040, 0.2790, 0.0170   ','X, Y, Z, SISSA 2020 + TP-AGB + VWPAGB')
    elif str(iop) == '-612':
        prihdr['modltrks'] = ('0.6960, 0.2840, 0.0200   ','X, Y, Z, SISSA 2020 + TP-AGB + VWPAGB')
    elif str(iop) == '-613':
        prihdr['modltrks'] = ('0.6680, 0.3020, 0.0300   ','X, Y, Z, SISSA 2020 + TP-AGB + VWPAGB')
    elif str(iop) == '-614':
        prihdr['modltrks'] = ('0.6390, 0.3210, 0.0400   ','X, Y, Z, SISSA 2020 + TP-AGB + VWPAGB')
    elif str(iop) == '-616':
        prihdr['modltrks'] = ('0.5840, 0.3560, 0.0600   ','X, Y, Z, SISSA 2020 + TP-AGB + VWPAGB')

    # CB2021 models, SISSA 2020 tracks + MBPAGB
    elif str(iop) == '-700':
        prihdr['modltrks'] = ('0.7700, 0.2300, 0.0000   ','X, Y, Z, SISSA 2020 + TP-AGB + MBPAGB')
    elif str(iop) == '-701':
        prihdr['modltrks'] = ('0.7509, 0.2490, 0.0001   ','X, Y, Z, SISSA 2020 + TP-AGB + MBPAGB')
    elif str(iop) == '-702':
        prihdr['modltrks'] = ('0.7508, 0.2490, 0.0002   ','X, Y, Z, SISSA 2020 + TP-AGB + MBPAGB')
    elif str(iop) == '-703':
        prihdr['modltrks'] = ('0.7505, 0.2490, 0.0005   ','X, Y, Z, SISSA 2020 + TP-AGB + MBPAGB')
    elif str(iop) == '-704':
        prihdr['modltrks'] = ('0.7490, 0.2500, 0.0010   ','X, Y, Z, SISSA 2020 + TP-AGB + MBPAGB')
    elif str(iop) == '-705':
        prihdr['modltrks'] = ('0.7460, 0.2520, 0.0020   ','X, Y, Z, SISSA 2020 + TP-AGB + MBPAGB')
    elif str(iop) == '-706':
        prihdr['modltrks'] = ('0.7400, 0.2560, 0.0040   ','X, Y, Z, SISSA 2020 + TP-AGB + MBPAGB')
    elif str(iop) == '-707':
        prihdr['modltrks'] = ('0.7350, 0.2590, 0.0060   ','X, Y, Z, SISSA 2020 + TP-AGB + MBPAGB')
    elif str(iop) == '-708':
        prihdr['modltrks'] = ('0.7290, 0.2630, 0.0080   ','X, Y, Z, SISSA 2020 + TP-AGB + MBPAGB')
    elif str(iop) == '-709':
        prihdr['modltrks'] = ('0.7230, 0.2670, 0.0100   ','X, Y, Z, SISSA 2020 + TP-AGB + MBPAGB')
    elif str(iop) == '-710':
        prihdr['modltrks'] = ('0.7130, 0.2730, 0.0140   ','X, Y, Z, SISSA 2020 + TP-AGB + MBPAGB')
    elif str(iop) == '-711':
        prihdr['modltrks'] = ('0.7040, 0.2790, 0.0170   ','X, Y, Z, SISSA 2020 + TP-AGB + MBPAGB')
    elif str(iop) == '-712':
        prihdr['modltrks'] = ('0.6960, 0.2840, 0.0200   ','X, Y, Z, SISSA 2020 + TP-AGB + MBPAGB')
    elif str(iop) == '-713':
        prihdr['modltrks'] = ('0.6680, 0.3020, 0.0300   ','X, Y, Z, SISSA 2020 + TP-AGB + MBPAGB')
    elif str(iop) == '-714':
        prihdr['modltrks'] = ('0.6390, 0.3210, 0.0400   ','X, Y, Z, SISSA 2020 + TP-AGB + MBPAGB')
    elif str(iop) == '-716':
        prihdr['modltrks'] = ('0.5840, 0.3560, 0.0600   ','X, Y, Z, SISSA 2020 + TP-AGB + MBPAGB')

    # BC2022 models, SISSA 2021 tracks + VWPAGB
    elif str(iop) == '-800':
        prihdr['modltrks'] = ('0.7700, 0.2300, 0.0000   ','X, Y, Z, SISSA 2021 + TP-AGB + VWPAGB')
    elif str(iop) == '-801':
        prihdr['modltrks'] = ('0.7509, 0.2490, 0.0001   ','X, Y, Z, SISSA 2021 + TP-AGB + VWPAGB')
    elif str(iop) == '-802':
        prihdr['modltrks'] = ('0.7508, 0.2490, 0.0002   ','X, Y, Z, SISSA 2021 + TP-AGB + VWPAGB')
    elif str(iop) == '-803':
        prihdr['modltrks'] = ('0.7505, 0.2490, 0.0005   ','X, Y, Z, SISSA 2021 + TP-AGB + VWPAGB')
    elif str(iop) == '-804':
        prihdr['modltrks'] = ('0.7490, 0.2500, 0.0010   ','X, Y, Z, SISSA 2021 + TP-AGB + VWPAGB')
    elif str(iop) == '-805':
        prihdr['modltrks'] = ('0.7460, 0.2520, 0.0020   ','X, Y, Z, SISSA 2021 + TP-AGB + VWPAGB')
    elif str(iop) == '-806':
        prihdr['modltrks'] = ('0.7400, 0.2560, 0.0040   ','X, Y, Z, SISSA 2021 + TP-AGB + VWPAGB')
    elif str(iop) == '-807':
        prihdr['modltrks'] = ('0.7350, 0.2590, 0.0060   ','X, Y, Z, SISSA 2021 + TP-AGB + VWPAGB')
    elif str(iop) == '-808':
        prihdr['modltrks'] = ('0.7290, 0.2630, 0.0080   ','X, Y, Z, SISSA 2021 + TP-AGB + VWPAGB')
    elif str(iop) == '-809':
        prihdr['modltrks'] = ('0.7230, 0.2670, 0.0100   ','X, Y, Z, SISSA 2021 + TP-AGB + VWPAGB')
    elif str(iop) == '-810':
        prihdr['modltrks'] = ('0.7130, 0.2730, 0.0140   ','X, Y, Z, SISSA 2021 + TP-AGB + VWPAGB')
    elif str(iop) == '-811':
        prihdr['modltrks'] = ('0.7040, 0.2790, 0.0170   ','X, Y, Z, SISSA 2021 + TP-AGB + VWPAGB')
    elif str(iop) == '-812':
        prihdr['modltrks'] = ('0.6960, 0.2840, 0.0200   ','X, Y, Z, SISSA 2021 + TP-AGB + VWPAGB')
    elif str(iop) == '-813':
        prihdr['modltrks'] = ('0.6680, 0.3020, 0.0300   ','X, Y, Z, SISSA 2021 + TP-AGB + VWPAGB')
    elif str(iop) == '-814':
        prihdr['modltrks'] = ('0.6390, 0.3210, 0.0400   ','X, Y, Z, SISSA 2021 + TP-AGB + VWPAGB')
    elif str(iop) == '-816':
        prihdr['modltrks'] = ('0.5840, 0.3560, 0.0600   ','X, Y, Z, SISSA 2021 + TP-AGB + VWPAGB')

    # CB2022 models, SISSA 2021 tracks + MBPAGB
    elif str(iop) == '-900':
        prihdr['modltrks'] = ('0.7700, 0.2300, 0.0000   ','X, Y, Z, SISSA 2021 + TP-AGB + MBPAGB')
    elif str(iop) == '-901':
        prihdr['modltrks'] = ('0.7509, 0.2490, 0.0001   ','X, Y, Z, SISSA 2021 + TP-AGB + MBPAGB')
    elif str(iop) == '-902':
        prihdr['modltrks'] = ('0.7508, 0.2490, 0.0002   ','X, Y, Z, SISSA 2021 + TP-AGB + MBPAGB')
    elif str(iop) == '-903':
        prihdr['modltrks'] = ('0.7505, 0.2490, 0.0005   ','X, Y, Z, SISSA 2021 + TP-AGB + MBPAGB')
    elif str(iop) == '-904':
        prihdr['modltrks'] = ('0.7490, 0.2500, 0.0010   ','X, Y, Z, SISSA 2021 + TP-AGB + MBPAGB')
    elif str(iop) == '-905':
        prihdr['modltrks'] = ('0.7460, 0.2520, 0.0020   ','X, Y, Z, SISSA 2021 + TP-AGB + MBPAGB')
    elif str(iop) == '-906':
        prihdr['modltrks'] = ('0.7400, 0.2560, 0.0040   ','X, Y, Z, SISSA 2021 + TP-AGB + MBPAGB')
    elif str(iop) == '-907':
        prihdr['modltrks'] = ('0.7350, 0.2590, 0.0060   ','X, Y, Z, SISSA 2021 + TP-AGB + MBPAGB')
    elif str(iop) == '-908':
        prihdr['modltrks'] = ('0.7290, 0.2630, 0.0080   ','X, Y, Z, SISSA 2021 + TP-AGB + MBPAGB')
    elif str(iop) == '-909':
        prihdr['modltrks'] = ('0.7230, 0.2670, 0.0100   ','X, Y, Z, SISSA 2021 + TP-AGB + MBPAGB')
    elif str(iop) == '-910':
        prihdr['modltrks'] = ('0.7130, 0.2730, 0.0140   ','X, Y, Z, SISSA 2021 + TP-AGB + MBPAGB')
    elif str(iop) == '-911':
        prihdr['modltrks'] = ('0.7040, 0.2790, 0.0170   ','X, Y, Z, SISSA 2021 + TP-AGB + MBPAGB')
    elif str(iop) == '-912':
        prihdr['modltrks'] = ('0.6960, 0.2840, 0.0200   ','X, Y, Z, SISSA 2021 + TP-AGB + MBPAGB')
    elif str(iop) == '-913':
        prihdr['modltrks'] = ('0.6680, 0.3020, 0.0300   ','X, Y, Z, SISSA 2021 + TP-AGB + MBPAGB')
    elif str(iop) == '-914':
        prihdr['modltrks'] = ('0.6390, 0.3210, 0.0400   ','X, Y, Z, SISSA 2021 + TP-AGB + MBPAGB')
    elif str(iop) == '-916':
        prihdr['modltrks'] = ('0.5840, 0.3560, 0.0600   ','X, Y, Z, SISSA 2021 + TP-AGB + MBPAGB')

    prihdr['modlimf']  = ('0.10, 100                ','[Mo], ML, MU, Chabrier (2003) IMF')
    prihdr['modltstp'] = ('1E4, 15E9, 221           ','[yr], Initial, Final age, Nsteps')
    prihdr['modlwavl'] = ('   5.6, 3.6E+8, var,13216','[A], Initial, Final wavl, Step, Npoints')
    prihdr['modlstl1'] = ('   5.6,  911.0, 0.9, 1007','Tlusty,Martins,UVBlue,WMBasic,PoWR,Rauch')
    prihdr['modlstl2'] = (' 911.5, 3540.5, 0.5, 5259','Tlusty,Martins,UVBlue,WMBasic,PoWR,Rauch')
    prihdr['modlstl3'] = ('3541.4, 7349.3, 0.9, 4232','Miles')
    prihdr['modlstl4'] = ('7351.0, 8750.0, 1.0, 1400','Stelib')
    prihdr['modlstl5'] = ('8770.0, 3.6E+8, var, 1318','BaSeL3.1, Aringer+, IRTF, Dusty TP-AGB')
    prihdr['modlbnt0']  =   ''
    prihdr['modlbnt1'] = ('BINTABLE 1                    ','SED of SSP vs. age in Lo/A')
    prihdr['modlbt1a'] = ('  SED(t)                      ','SED for 221 values of age t in Lo/A')
    prihdr['modlbt1b'] = ('  SED(t=0) = SED(t=10 Myr)    ','First time step, n=1')
    prihdr['modlbt1c'] = ('  SED(t=10 Myr)               ','Second time step, n=2')
    prihdr['modlbt1d'] = ('  SED(t=15 Gyr)               ','Last time step, n=221')
    prihdr['modlbt20']  =   ''
    prihdr['modlbnt2'] = ('BINTABLE 2                    ','PHYSICAL PROPERTIES OF SSP vs. AGE')
    prihdr['modlbt2a'] = ('  logNLy, logNHeI, logNHeII   ','H, HeI, and HeII ionizing photons')
    prihdr['modlbt2b'] = ('  B912                        ','Amplitude of Lyman break')
    prihdr['modlbt2c'] = ('  B4000VN, B4000SDSS, B4000   ','Amplitude of 4000 A break')
    prihdr['modlbt2d'] = ('  SNR_yr                      ','Super nova rate /yr')
    prihdr['modlbt2e'] = ('  PNBR_yr                     ','Planetary Nebula birth rate /yr')
    prihdr['modlbt2f'] = ('  NBH                         ','Number of Black Holes formed')
    prihdr['modlbt2g'] = ('  NNS                         ','Number of Neutron Stars formed')
    prihdr['modlbt2h'] = ('  NWD                         ','Number of White Dwarfs formed')
    prihdr['modlbt2i'] = ('  Lx(Lo)                      ','X ray luminosity (Lo)')
    prihdr['modlbt2j'] = ('  Mstars                      ','Mass in living stars in Mo')
    prihdr['modlbt2k'] = ('  Mremnants                   ','Mass in stellar remnants in Mo')
    prihdr['modlbt2l'] = ('  Mretgas                     ','Mass of gas returned to ISM in Mo')
    prihdr['modlbt2m'] = ('  Mgalaxy                     ','Mass of galaxy Mo')
    prihdr['modlbt2n'] = ('  SFR_yr                      ','Star Formation Rate in Mo/yr')
    prihdr['modlbt2o'] = ('  Mtot                        ','Mstars + Mremnants in Mo')
    prihdr['modlbt2p'] = ('  Mtot_Lb, Mtot_Lv, Mtot_Lk   ','Mtot/light ratio in B,V,K bands')
    prihdr['modlbt2q'] = ('  Mliv_Lb, Mliv_Lv, Mliv_Lk   ','Mstars/light ratio in B,V,K bands')
    prihdr['modlbt2r'] = ('  bt_yr                       ','Evolutionary flux,Renzini & Buzzoni')
    prihdr['modlbt2s'] = ('  Bt_yr_Lo                    ','b(t)/Bolometric luminosity.')
    prihdr['modlbt2t'] = ('  TurnoffMass                 ','Turnoff Mass in Mo')
    prihdr['modlbt2u'] = ('  BPMS_BMS                    ','Ratio of PostMS/MS bolometric flux')
    prihdr['modlbt30']  =   ''
    prihdr['modlbnt3'] = ('BINTABLE 3                    ','MAGNITUDES vs. AGE')
    prihdr['modlbt3a'] = ('  VEGA mags:                  ','Vega magnitude system')
    prihdr['modlbt3b'] = ('     Mbol:                    ','Bolometric magnitude')
    prihdr['modlbt3c'] = ('     UBVRIJKL:                ','Johnson filters')
    prihdr['modlbt3d'] = ('     RI:                      ','Cousins filters')
    prihdr['modlbt3e'] = ('     JHK:                     ','Palomar filters')
    prihdr['modlbt3f'] = ('     Kprime:                  ','Cowie Kprime filter')
    prihdr['modlbt3g'] = ('     JHKs:                    ','2Mass filters')
    prihdr['modlbt3h'] = ('     3.6,4.5,5.7,7.9 micron:  ','IRAC filters')
    prihdr['modlbt3i'] = ('     12,25,60,100 micron:     ','IRAS filters')
    prihdr['modlbt3j'] = ('     24,70,160 micron:        ','MIPS filters')
    prihdr['modlbt3k'] = ('     F220w,F250w,F330w,F410w, ','HST ACS WFC wide filters')
    prihdr['modlbt3l'] = ('     F435w,F475w,F555w,F606w, ',' "   "   "   "      "   ')
    prihdr['modlbt3m'] = ('     F625w,F775w,F814w:       ',' "   "   "   "      "   ')
    prihdr['modlbt3n'] = ('     f225w,f275w,f336w,f438w, ','HST UVIS1 filters')
    prihdr['modlbt3o'] = ('     f547m,f555w,f606w,f625w, ',' "    "      "   ')
    prihdr['modlbt3p'] = ('     f656n,f657n,f658n,f814w: ',' "    "      "   ')
    prihdr['modlbt3q'] = ('  AB mags:                    ','AB magnitude system')
    prihdr['modlbt3r'] = ('     ugriz:                   ','SDSS filters')
    prihdr['modlbt3s'] = ('     ugriyzKs:                ','CFHT MegaCam filters')
    prihdr['modlbt3t'] = ('     FUV,NUV:                 ','GALEX filters')
    prihdr['modlbt40']  =   ''
    prihdr['modlbnt4'] = ('BINTABLE 4                    ','LINE STRENGTH INDICES vs. AGE')
    prihdr['modlbt4a'] = ('  CN1,CN2,Ca4227,G4300,Fe4383,','Lick Indices')
    prihdr['modlbt4b'] = ('  Ca4455,Fe4531,Fe4668,Hbeta, ','  "     "')
    prihdr['modlbt4c'] = ('  Fe5015,Mg1,Mg2,Mgb,Fe5270,  ','  "     "')
    prihdr['modlbt4d'] = ('  Fe5335,Fe5406,Fe5709,Fe5782,','  "     "')
    prihdr['modlbt4e'] = ('  NaD,TiO1,TiO2               ','  "     "')
    prihdr['modlbt4f'] = ('  HdelA, HgamA, HdelF, HgamF  ','Worthey & Ottaviani (1997)')
    prihdr['modlbt4g'] = ('  CaII8498, CaII8542, CaII8662','Diaz, Terlevich, &')
    prihdr['modlbt4h'] = ('  MgI8807                     ','      Terlevich (1989)')
    prihdr['modlbt4i'] = ('  H8_3889, H9_3835, H10_3798  ','Marcillac et al. (2006)')
    prihdr['modlbt4j'] = ('  BHHK                        ','Brodie & Hanes HK index')
    prihdr['modlbt4k'] = ('  D4000                       ','4000 A break index (Gorgas&Cardiel)')
    prihdr['modlbt4l'] = ('  B4000VN                     ','Very narrow 4000 A break')
    prihdr['modlbt4m'] = ('  BL1302,SiIV,BL1425,Fe1453,  ','UV, Fanelli et al. (1987,1990,1992)')
    prihdr['modlbt4n'] = ('  CIV1548a,CIV1548c,CIV1548e, ',' "     "      "         "')
    prihdr['modlbt4o'] = ('  BL1617,BL1664,BL1719,BL1853,',' "     "      "         "')
    prihdr['modlbt4p'] = ('  FeII2402,BL2538,FeII2609,   ',' "     "      "         "')
    prihdr['modlbt4q'] = ('  MgII,MgI,Mgwide,FeI,BL3096  ',' "     "      "         "')
    prihdr['modlbt4r'] = ('  CIVa                        ','CIV  1540 A stellar absorption')
    prihdr['modlbt4s'] = ('  HeIIe                       ','HeII 1640 A stellar emission')

    # Build fits compatible header for BINTABLE2 (physical properties)
    sechdr = hdul[2].header
    sechdr['TTYPE1']  = 'logage'             # label for column 1
    sechdr['TTYPE2']  = 'logNLy'             # label for column 2
    sechdr['TTYPE3']  = 'logNHeI'            # label for column 3
    sechdr['TTYPE4']  = 'logNHeII'           # label for column 4
    sechdr['TTYPE5']  = 'logNHeII-logNLy'    # label for column 5
    sechdr['TTYPE6']  = 'B912'               # label for column 6
    sechdr['TTYPE7']  = 'B4000VN'            # label for column 7
    sechdr['TTYPE8']  = 'B4000SDSS'          # label for column 8
    sechdr['TTYPE9']  = 'B4000'              # label for column 9

    sechdr['TTYPE10'] = 'SNR'                # label for column 10
    sechdr['TTYPE11'] = 'PISNR'              # label for column 11
    sechdr['TTYPE12'] = 'RegSNR'             # label for column 12
    sechdr['TTYPE13'] = 'TypeIaSNR'          # label for column 13
    sechdr['TTYPE14'] = 'FailedSNR'          # label for column 14
    sechdr['TTYPE15'] = 'PNBR'               # label for column 15
    sechdr['TTYPE16'] = 'NBH'                # label for column 16
    sechdr['TTYPE17'] = 'NNS'                # label for column 17
    sechdr['TTYPE18'] = 'NWD'                # label for column 18
    sechdr['TTYPE19'] = 'Lx'                 # label for column 19

    sechdr['TTYPE20'] = 'Mstars'             # label for column 20
    sechdr['TTYPE21'] = 'Mremnants'          # label for column 21
    sechdr['TTYPE22'] = 'Mretgas'            # label for column 22
    sechdr['TTYPE23'] = 'Mgalaxy'            # label for column 23
    sechdr['TTYPE24'] = 'SFR'                # label for column 24
    sechdr['TTYPE25'] = 'Mtot'               # label for column 25
    sechdr['TTYPE26'] = 'Mtot_Lb'            # label for column 26
    sechdr['TTYPE27'] = 'Mtot_Lv'            # label for column 27
    sechdr['TTYPE28'] = 'Mtot_Lk'            # label for column 28
    sechdr['TTYPE29'] = 'Mliv_Lb'            # label for column 29

    sechdr['TTYPE30'] = 'Mliv_Lv'            # label for column 30
    sechdr['TTYPE31'] = 'Mliv_Lk'            # label for column 31
    sechdr['TTYPE32'] = 'EvolutionaryFlux'   # label for column 32
    sechdr['TTYPE33'] = 'SpecificEvolFlux'   # label for column 33
    sechdr['TTYPE34'] = 'TurnoffMass'        # label for column 34
    sechdr['TTYPE35'] = 'BPMS_BMS'           # label for column 35
    sechdr['TTYPE36'] = 'TotalMassLossRate'  # label for column 36
    sechdr['TTYPE37'] = 'DustProductionRate' # label for column 37

    sechdr['TUNIT1']  = 'yr            (log age of SSP)'
    sechdr['TUNIT2']  = 'photons/Mo    (log number of HI ionizing photons)'
    sechdr['TUNIT3']  = 'photons/Mo    (log number of HeI ionizing photons)'
    sechdr['TUNIT4']  = 'photons/Mo    (log number of HeII ionizing photons)'
    sechdr['TUNIT5']  = '..........    (log ratio number of ionizing HeII/HI)'
    sechdr['TUNIT6']  = '..........    (Lyman break amplitude)'
    sechdr['TUNIT7']  = '..........    (4000 A break, very narrow definition)'
    sechdr['TUNIT8']  = '..........    (4000 A break, SDSS definition)'
    sechdr['TUNIT9']  = '..........    (4000 A break, original definition)'

    sechdr['TUNIT10'] = 'SN/yr         (Total Supernova rate)'
    sechdr['TUNIT11'] = 'SN/yr         (PI Supernova rate)'
    sechdr['TUNIT12'] = 'SN/yr         (Regular Supernova rate)'
    sechdr['TUNIT13'] = 'SN/yr         (Type Ia Supernova rate)'
    sechdr['TUNIT14'] = 'SN/yr         (Failed Supernova rate)'
    sechdr['TUNIT15'] = 'PN/yr         (Planetary Nebula birth rate)'
    sechdr['TUNIT16'] = 'number/Mo     (Number of Black holes)'
    sechdr['TUNIT17'] = 'number/Mo     (Number of Neutron stars)'
    sechdr['TUNIT18'] = 'number/Mo     (Number of White dwarfs)'
    sechdr['TUNIT19'] = 'Lo            (X ray luminosity)'

    sechdr['TUNIT20'] = 'Mo            (Mliv = Mass in living stars)'
    sechdr['TUNIT21'] = 'Mo            (Mrem = Mass in stellar remnants)'
    sechdr['TUNIT22'] = 'Mo            (Mgas = Mass of processed gas)'
    sechdr['TUNIT23'] = 'Mo            (Mass of galaxy)'
    sechdr['TUNIT24'] = 'Mo/yr         (Star formation rate)'
    sechdr['TUNIT25'] = 'Mo            (Mtot = Mliv + Mrem)'
    sechdr['TUNIT26'] = 'Mo/Lo         (Mtot/Light ratio B band)'
    sechdr['TUNIT27'] = 'Mo/Lo         (Mtot/Light ratio V band)'
    sechdr['TUNIT28'] = 'Mo/Lo         (Mtot/Light ratio K band)'
    sechdr['TUNIT29'] = 'Mo/Lo         (Mliv/Light ratio B band)'

    sechdr['TUNIT30'] = 'Mo/Lo         (Mliv/Light ratio V band)'
    sechdr['TUNIT31'] = 'Mo/Lo         (Mliv/Light ratio K band)'
    sechdr['TUNIT32'] = 'Nstars/yr     (Evolutionary flux)'
    sechdr['TUNIT33'] = 'Nstars/yr/Lo  (Specific Evolutionary flux)'
    sechdr['TUNIT34'] = 'Mo            (MS Turnoff mass)'
    sechdr['TUNIT35'] = '..........    (PostMS/MS bolometric flux)'
    sechdr['TUNIT36'] = 'M/Mo/yr       (Total Mass Loss Rate)'
    sechdr['TUNIT37'] = 'M/Mo/yr       (Total Dust Production Rate)'

    # Build fits compatible header for BINTABLE3 (photometry)
    sechdr = hdul[3].header
    sechdr['TTYPE1']  = 'logage'
    sechdr['TUNIT1']  = 'yr       (log age of SED)'        # label for column 1
    sechdr['TUNIT2']  = 'Vega mag (Bolometric magnitude)'  # label for column 2
    sechdr['TUNIT3']  = 'AB mag   (FUV_GALEX_AB)'          # label for column 3
    sechdr['TUNIT4']  = 'AB mag   (ACSWFC_F220w_AB)'       # label for column 4
    sechdr['TUNIT5']  = 'AB mag   (NUV_GALEX_AB)'          # label for column 5
    sechdr['TUNIT6']  = 'AB mag   (WFC3_F225W_AB)'         # label for column 6
    sechdr['TUNIT7']  = 'AB mag   (UVIS1_f225w_AB)'        # label for column 7
    sechdr['TUNIT8']  = 'AB mag   (UVIS1_f275w_AB)'        # label for column 8
    sechdr['TUNIT9']  = 'AB mag   (ACSWFC_F250w_AB)'       # label for column 9
    sechdr['TUNIT10'] = 'AB mag   (ACSWFC_F330w_AB)'       # label for column 10
    sechdr['TUNIT11'] = 'AB mag   (UVIS1_f336w_AB)'        # label for column 11
    sechdr['TUNIT12'] = 'AB mag   (WFC3_F336W_AB)'         # label for column 12
    sechdr['TUNIT13'] = 'AB mag   (u_SDSS_AB)'             # label for column 13
    sechdr['TUNIT14'] = 'Vega mag (U_Johnson)'             # label for column 14
    sechdr['TUNIT15'] = 'AB mag   (u3_CFHT_MC_AB)'         # label for column 15
    sechdr['TUNIT16'] = 'AB mag   (u1_CFHT_MC_AB)'         # label for column 16
    sechdr['TUNIT17'] = 'AB mag   (ACSWFC_F410w_AB)'       # label for column 17
    sechdr['TUNIT18'] = 'AB mag   (WFC3_FR388N_AB)'        # label for column 18
    sechdr['TUNIT19'] = 'AB mag   (ACSWFC_F435w_AB)'       # label for column 19
    sechdr['TUNIT20'] = 'AB mag   (WFC3_F438W_AB)'         # label for column 20
    sechdr['TUNIT21'] = 'AB mag   (UVIS1_f438w_AB)'        # label for column 21
    sechdr['TUNIT22'] = 'Vega mag (B3_Johnson)'            # label for column 22
    sechdr['TUNIT23'] = 'Vega mag (B2_Johnson)'            # label for column 23
    sechdr['TUNIT24'] = 'AB mag   (g_SDSS_AB)'             # label for column 24
    sechdr['TUNIT25'] = 'AB mag   (ACSWFC_F475w_AB)'       # label for column 25
    sechdr['TUNIT26'] = 'AB mag   (g3_CFHT_MC_AB)'         # label for column 26
    sechdr['TUNIT27'] = 'AB mag   (g1_CFHT_MC_AB)'         # label for column 27
    sechdr['TUNIT28'] = 'AB mag   (UVIS1_f555w_AB)'        # label for column 28
    sechdr['TUNIT29'] = 'AB mag   (WFC3_F555W_AB)'         # label for column 29
    sechdr['TUNIT30'] = 'AB mag   (ACSWFC_F555w_AB)'       # label for column 30
    sechdr['TUNIT31'] = 'AB mag   (UVIS1_f547m_AB)'        # label for column 31
    sechdr['TUNIT32'] = 'Vega mag (V_Johnson)'             # label for column 32
    sechdr['TUNIT33'] = 'AB mag   (ACSWFC_F606w_AB)'       # label for column 33
    sechdr['TUNIT34'] = 'AB mag   (UVIS1_f606w_AB)'        # label for column 34
    sechdr['TUNIT35'] = 'AB mag   (r_SDSS_AB)'             # label for column 35
    sechdr['TUNIT36'] = 'AB mag   (r1_CFHT_MC_AB)'         # label for column 36
    sechdr['TUNIT37'] = 'AB mag   (UVIS1_f625w_AB)'        # label for column 37
    sechdr['TUNIT38'] = 'AB mag   (ACSWFC_F625w_AB)'       # label for column 38
    sechdr['TUNIT39'] = 'AB mag   (r3_CFHT_MC_AB)'         # label for column 39
    sechdr['TUNIT40'] = 'AB mag   (UVIS1_f656n_AB)'        # label for column 40
    sechdr['TUNIT41'] = 'AB mag   (UVIS1_f657n_AB)'        # label for column 41
    sechdr['TUNIT42'] = 'AB mag   (UVIS1_f658n_AB)'        # label for column 42
    sechdr['TUNIT43'] = 'Vega mag (R_Cousins)'             # label for column 43
    sechdr['TUNIT44'] = 'Vega mag (R_Johnson)'             # label for column 44
    sechdr['TUNIT45'] = 'AB mag   (i_SDSS_AB)'             # label for column 45
    sechdr['TUNIT46'] = 'AB mag   (i2_CFHT_MC_AB)'         # label for column 46
    sechdr['TUNIT47'] = 'AB mag   (i3_CFHT_MC_AB)'         # label for column 47
    sechdr['TUNIT48'] = 'AB mag   (ACSWFC_F775w_AB)'       # label for column 48
    sechdr['TUNIT49'] = 'Vega mag (I_Cousins)'             # label for column 49
    sechdr['TUNIT50'] = 'AB mag   (UVIS1_f814w_AB)'        # label for column 50
    sechdr['TUNIT51'] = 'AB mag   (ACSWFC_F814w_AB)'       # label for column 51
    sechdr['TUNIT52'] = 'AB mag   (WFC3_F814W_AB)'         # label for column 52
    sechdr['TUNIT53'] = 'Vega mag (I_Johnson)'             # label for column 53
    sechdr['TUNIT54'] = 'AB mag   (z1_CFHT_MC_AB)'         # label for column 54
    sechdr['TUNIT55'] = 'AB mag   (z_SDSS_AB)'             # label for column 55
    sechdr['TUNIT56'] = 'AB mag   (z3_CFHT_MC_AB)'         # label for column 56
    sechdr['TUNIT57'] = 'AB mag   (WFC3_F110W_AB)'         # label for column 57
    sechdr['TUNIT58'] = 'Vega mag (J_2Mass)'               # label for column 58
    sechdr['TUNIT59'] = 'Vega mag (J_Johnson)'             # label for column 59
    sechdr['TUNIT60'] = 'AB mag   (WFC3_F125W_AB)'         # label for column 60
    sechdr['TUNIT61'] = 'Vega mag (J_Palomar)'             # label for column 61
    sechdr['TUNIT62'] = 'AB mag   (WFC3_F160W_AB)'         # label for column 62
    sechdr['TUNIT63'] = 'Vega mag (H_Palomar)'             # label for column 63
    sechdr['TUNIT64'] = 'AB mag   (H_2MASS_AB)'            # label for column 64
    sechdr['TUNIT65'] = 'Vega mag (H_2Massr)'              # label for column 65
    sechdr['TUNIT66'] = 'Vega mag (K    prime_Cowie)'          # label for column 66
    sechdr['TUNIT67'] = 'AB mag   (Ks_CFHT_WC_AB)'         # label for column 67
    sechdr['TUNIT68'] = 'Vega mag (Ks_2Mass)'              # label for column 68
    sechdr['TUNIT69'] = 'Vega mag (K_Johnson)'             # label for column 69
    sechdr['TUNIT70'] = 'Vega mag (K_Palomar)'             # label for column 70
    sechdr['TUNIT71'] = 'Vega mag (L_Johnson)'             # label for column 71
    sechdr['TUNIT72'] = 'Vega mag (I3p6_IRAC)'             # label for column 72
    sechdr['TUNIT73'] = 'Vega mag (I4p5_IRAC)'             # label for column 73
    sechdr['TUNIT74'] = 'Vega mag (I5p7_IRAC)'             # label for column 74
    sechdr['TUNIT75'] = 'Vega mag (I7p9_IRAC)'             # label for column 75
    sechdr['TUNIT76'] = 'Vega mag (I12_IRAS)'              # label for column 76
    sechdr['TUNIT77'] = 'Vega mag (M24_MIPS)'              # label for column 77
    sechdr['TUNIT78'] = 'Vega mag (I25_IRAS)'              # label for column 78
    sechdr['TUNIT79'] = 'Vega mag (I60_IRAS)'              # label for column 79
    sechdr['TUNIT80'] = 'Vega mag (M70_MIPS)'              # label for column 80
    sechdr['TUNIT81'] = 'Vega mag (I100_IRAS)'             # label for column 81
    sechdr['TUNIT82'] = 'Vega mag (M160_MIPS)'             # label for column 82
    sechdr['TUNIT83'] = 'FUV Flux (GALEX)'                 # label for column 83
    sechdr['TUNIT84'] = 'NUV Flux (GALEX)'                 # label for column 84
    sechdr['TUNIT85'] = '1500 A Flux (square)'             # label for column 85

    # Build fits compatible header for BINTABLE4 (line indices)
    sechdr = hdul[4].header
    sechdr['TTYPE1']  = 'logage'             # label for column 1
    sechdr['TTYPE2']  = 'CN1'                # label for column 2
    sechdr['TTYPE3']  = 'CN2'                # label for column 3
    sechdr['TTYPE4']  = 'Ca4227'             # label for column 4
    sechdr['TTYPE5']  = 'G4300'              # label for column 5
    sechdr['TTYPE6']  = 'Fe4383'             # label for column 6
    sechdr['TTYPE7']  = 'Ca4455'             # label for column 7
    sechdr['TTYPE8']  = 'Fe4531'             # label for column 8
    sechdr['TTYPE9']  = 'Fe4668'             # label for column 9
    sechdr['TTYPE10'] = 'Hbeta'              # label for column 10
    sechdr['TTYPE11'] = 'Fe5015'             # label for column 11
    sechdr['TTYPE12'] = 'Mg1'                # label for column 12
    sechdr['TTYPE13'] = 'Mg2'                # label for column 13
    sechdr['TTYPE14'] = 'Mgb'                # label for column 14
    sechdr['TTYPE15'] = 'Fe5270'             # label for column 15
    sechdr['TTYPE16'] = 'Fe5335'             # label for column 16
    sechdr['TTYPE17'] = 'Fe5406'             # label for column 17
    sechdr['TTYPE18'] = 'Fe5709'             # label for column 18
    sechdr['TTYPE19'] = 'Fe5782'             # label for column 19
    sechdr['TTYPE20'] = 'NaD'                # label for column 20
    sechdr['TTYPE21'] = 'TiO1'               # label for column 21
    sechdr['TTYPE22'] = 'TiO2'               # label for column 22
    sechdr['TTYPE23'] = 'HdeltaA'            # label for column 23
    sechdr['TTYPE24'] = 'HgammaA'            # label for column 24
    sechdr['TTYPE25'] = 'HdeltaF'            # label for column 25
    sechdr['TTYPE26'] = 'HgammaF'            # label for column 26
    sechdr['TTYPE27'] = 'D4000'              # label for column 27
    sechdr['TTYPE28'] = 'B4000VN'            # label for column 28
    sechdr['TTYPE29'] = 'CaII8498'           # label for column 29
    sechdr['TTYPE30'] = 'CaII8542'           # label for column 30
    sechdr['TTYPE31'] = 'CaII8662'           # label for column 31
    sechdr['TTYPE32'] = 'MgI8807'            # label for column 32
    sechdr['TTYPE33'] = 'H83889'             # label for column 33
    sechdr['TTYPE34'] = 'H93835'             # label for column 34
    sechdr['TTYPE35'] = 'H103798'            # label for column 35
    sechdr['TTYPE36'] = 'BHHK'               # label for column 36
    sechdr['TTYPE37'] = 'BL1302'             # label for column 37
    sechdr['TTYPE38'] = 'SiIV'               # label for column 38
    sechdr['TTYPE39'] = 'BL1425'             # label for column 39
    sechdr['TTYPE40'] = 'Fe1453'             # label for column 40
    sechdr['TTYPE41'] = 'CIV1548a'           # label for column 41
    sechdr['TTYPE42'] = 'CIV1548c'           # label for column 42
    sechdr['TTYPE43'] = 'CIV1548e'           # label for column 43
    sechdr['TTYPE44'] = 'BL1617'             # label for column 44
    sechdr['TTYPE45'] = 'BL1664'             # label for column 45
    sechdr['TTYPE46'] = 'BL1719'             # label for column 46
    sechdr['TTYPE47'] = 'BL1853'             # label for column 47
    sechdr['TTYPE48'] = 'FeII2402'           # label for column 48
    sechdr['TTYPE49'] = 'BL2538'             # label for column 49
    sechdr['TTYPE50'] = 'FeII2609'           # label for column 50
    sechdr['TTYPE51'] = 'MgII'               # label for column 51
    sechdr['TTYPE52'] = 'MgI'                # label for column 52
    sechdr['TTYPE53'] = 'Mgwide'             # label for column 53
    sechdr['TTYPE54'] = 'FeI'                # label for column 54
    sechdr['TTYPE55'] = 'BL3096'             # label for column 55
    sechdr['TTYPE56'] = 'CIVabs'             # label for column 56
    sechdr['TTYPE57'] = 'HeIIems'            # label for column 57
    sechdr['TUNIT1']  = 'yr   (log age of SED)'
    sechdr['TUNIT2']  = 'mag  (CN1 Lick index)'
    sechdr['TUNIT3']  = 'mag  (CN2 Lick index)'
    sechdr['TUNIT4']  = 'A    (Ca4227 Lick index)'
    sechdr['TUNIT5']  = 'A    (G4300 Lick index)'
    sechdr['TUNIT6']  = 'A    (Fe4383 Lick index)'
    sechdr['TUNIT7']  = 'A    (Ca4455 Lick index)'
    sechdr['TUNIT8']  = 'A    (Fe4531 Lick index)'
    sechdr['TUNIT9']  = 'A    (Fe4668 Lick index)'
    sechdr['TUNIT10'] = 'A    (Hbeta Lick index)'
    sechdr['TUNIT11'] = 'A    (Fe5015 Lick index)'
    sechdr['TUNIT12'] = 'mag  (Mg1 Lick index)'
    sechdr['TUNIT13'] = 'mag  (Mg2 Lick index)'
    sechdr['TUNIT14'] = 'A    (Mgb Lick index)'
    sechdr['TUNIT15'] = 'A    (Fe5270 Lick index)'
    sechdr['TUNIT16'] = 'A    (Fe5335 Lick index)'
    sechdr['TUNIT17'] = 'A    (Fe5406 Lick index)'
    sechdr['TUNIT18'] = 'A    (Fe5709 Lick index)'
    sechdr['TUNIT19'] = 'A    (Fe5782 Lick index)'
    sechdr['TUNIT20'] = 'A    (NaD Lick index)'
    sechdr['TUNIT21'] = 'mag  (TiO1 Lick index)'
    sechdr['TUNIT22'] = 'mag  (TiO2 Lick index)'
    sechdr['TUNIT23'] = 'A    (Worthey & Ottaviani HdeltaA index)'
    sechdr['TUNIT24'] = 'A    (Worthey & Ottaviani HgammaA index)'
    sechdr['TUNIT25'] = 'A    (Worthey & Ottaviani HdeltaF index)'
    sechdr['TUNIT26'] = 'A    (Worthey & Ottaviani HgammaF index)'
    sechdr['TUNIT27'] = '     (Gorgas & Cardiel 4000 A break index)'
    sechdr['TUNIT28'] = '     (Very narrow 4000 A break index)'
    sechdr['TUNIT29'] = 'A    (Diaz & Terlevich**2 CaII8498 index)'
    sechdr['TUNIT30'] = 'A    (Diaz & Terlevich**2 CaII8542 index)'
    sechdr['TUNIT31'] = 'A    (Diaz & Terlevich**2 CaII8662 index)'
    sechdr['TUNIT32'] = 'A    (Diaz & Terlevich**2 MgI8807 index)'
    sechdr['TUNIT33'] = 'A    (Marcillac et al. H83889 index)'
    sechdr['TUNIT34'] = 'A    (Marcillac et al. H93835 index)'
    sechdr['TUNIT35'] = 'A    (Marcillac et al. H103798 index)'
    sechdr['TUNIT36'] = 'A    (Brodie & Hanes HK index)'
    sechdr['TUNIT37'] = 'A    (Fanelli et al. BL1302 index)'
    sechdr['TUNIT38'] = 'A    (Fanelli et al. SiIV index)'
    sechdr['TUNIT39'] = 'A    (Fanelli et al. BL1425 index)'
    sechdr['TUNIT40'] = 'A    (Fanelli et al. Fe1453 index)'
    sechdr['TUNIT41'] = 'A    (Fanelli et al. CIV1548a index)'
    sechdr['TUNIT42'] = 'A    (Fanelli et al. CIV1548c index)'
    sechdr['TUNIT43'] = 'A    (Fanelli et al. CIV1548e index)'
    sechdr['TUNIT44'] = 'A    (Fanelli et al. BL1617 index)'
    sechdr['TUNIT45'] = 'A    (Fanelli et al. BL1664 index)'
    sechdr['TUNIT46'] = 'A    (Fanelli et al. BL1719 index)'
    sechdr['TUNIT47'] = 'A    (Fanelli et al. BL1853 index)'
    sechdr['TUNIT48'] = 'A    (Fanelli et al. FeII2402 index)'
    sechdr['TUNIT49'] = 'A    (Fanelli et al. BL2538 index)'
    sechdr['TUNIT50'] = 'A    (Fanelli et al. FeII2609 index)'
    sechdr['TUNIT51'] = 'A    (Fanelli et al. MgII index)'
    sechdr['TUNIT52'] = 'A    (Fanelli et al. MgI index)'
    sechdr['TUNIT53'] = 'A    (Fanelli et al. Mgwide index)'
    sechdr['TUNIT54'] = 'A    (Fanelli et al. FeI index)'
    sechdr['TUNIT55'] = 'A    (Fanelli et al. BL3096 index)'
    sechdr['TUNIT56'] = 'A    (CIV 1540 A stellar absorption line index)'
    sechdr['TUNIT57'] = 'A    (HeII 1640 A stellar emission line index)'

    def updhdr():
        from datetime import datetime
        global hdr
        # Updates model file header according to chosen SFR
        # ('#',10x,'I      S.F.R.: SSP = Zero Length Burst at t = 0',74X,'I')
        # ('#',10x,'I      S.F.R.: Exponential with MU9 = ',F5.3,'/Gyr',6X,'TAU = ',F7.3,' Gyr (includes processed gas recycling)',16X,'I')
        # ('#',10x,'I      S.F.R.: Finite Burst of Duration = ',1pe9.3,' yr',67X,'I')
        # ('#',10x,'I      S.F.R.: Constant = ',1pe10.3,' Mo/yr',79X,'I')
        # ('#',10x,'I      S.F.R.: Exponential with MU9 = ',F5.3,'/Gyr',6X,'TAU = ',F7.3,' Gyr (does not include processed gas recycling)', 8X,'I')
        # ('#',10x,'I      S.F.R.: 2 Bursts: Burst 1 at ',F6.3,'/Gyr, Strength = ',F8.5,', in file = ',a,10x,'I')
        # ('#',10x,'I                        Burst 2 at ',F6.3,'/Gyr, Strength = ',F8.5,', in file = ',a,10x,'I')
        # ('#',10x,'I      S.F.R.: Double Exponential with TAU_1 = ',F7.3,' Gyr (after Chen et al., does not include processed gas recycling) I')
        # ('#',10x,'I                                      TAU_2 = ',F7.3,' Gyr                                                               I')
        # ('#',10x,'I               which join smoothly at TAU_J = ',f7.3,' Gyr                                                               I')
        # ('#',10x,'I                    Look-back-time = Tform  = ',f7.3,' Gyr                                                               I')
        # ('#',10x,'I                            Burst amplitude = ',f7.3,'     ( = ratio of mass in burst / subyacent mass at t = Tform)     I')
        # ('#',10x,'I                            Burst starts at = ',f7.3,' Gyr                                                               I')
        # ('#',10x,'I                             Burst duration = ',f7.3,' Gyr                                                               I')

        # Update SFR info
        v=' '
        h2a = '' ; h2b = ''
        io  = bt.io
        ff  = bs.lnam(name)
        if io == 0:
            h1 = 'I'+6*v+'S.F.R.: SSP = Zero Length Burst at t = 0'
        elif io == 1:
            v1 = bt.hdr[0] ; v2 = bt.hdr[1]
            h1 = 'I'+6*v+'S.F.R.: Exponential with MU9 = '+f'{v1:5.3f}/Gyr'+6*v+'TAU = '+f'{v2*1e-9:7.3f} Gyr (includes processed gas recycling)'
        elif io == 2:
            v1 = bt.hdr[0]
            h1 = 'I'+6*v+'S.F.R.: Finite Burst of Duration = '+f'{v1:9.3E} yr'
        elif io == 3:
            v1 = bt.hdr[0]
            h1 = 'I'+6*v+'S.F.R.: Constant = '+f'{v1:10.3E} Mo/yr'
        elif io == 4:
            v1 = bt.hdr[0] ; v2 = bt.hdr[1]
            h1 = 'I'+6*v+'S.F.R.: Exponential with MU9 = '+f'{v1:5.3f}/Gyr'+6*v+'TAU = '+f'{v2*1e-9:7.3f} Gyr (does not include processed gas recycling)'
        elif io == 5:
            bp1 = bt.bp1 ; bp2 = bt.bp2
            v1 = bp1[0]*1.E-9 ; v2 = bp1[1] ; v3 = bs.lnam(bp1[2])
            h1 = 'I'+6*v+'S.F.R.: 2 Bursts: Burst 1 at ' + f'{v1:6.3f} Gyr, Strength = {v2:8.4f}, in file = ' + v3
            v1 = bp2[0]*1.E-9 ; v2 = bp2[1] ; v3 = bs.lnam(bp2[2])
            h2 = 'I'+6*v+'                  Burst 2 at ' + f'{v1:6.3f} Gyr, Strength = {v2:8.4f}, in file = ' + v3
            h2  = h2 + (121 - len(h2))*v + 'I'
            h2  = h2.replace(' ','|')
            h2a = h2[:61]
            h2b = h2[61:]
        elif io == 6:
            v1 = bt.hdr[0]
            h1 = 'I'+6*v+'S.F.R.: Delayed with TAU = '+f'{v1*1e-9:7.3f} Gyr'
        elif io == 7:
            h1 = 'I'+6*v+'S.F.R.: Read from table in file ' + bt.hdr[0]
        elif io == 8:
            v1 = bt.hdr[0]
            h1 = 'I'+6*v+'S.F.R.: Linear with TAU = '+f'{v1*1e-9:7.3f} Gyr'
        elif io == 9:
            h1 = 'I'+6*v+'S.F.R.: Follows Chen et al. (2012, MNRAS, 421, 314) with parameters in file ' + ff.replace('fits','sfr')
        h1  = h1 + (121 - len(h1))*v + 'I'
        h1  = h1.replace(' ','|')
        h1a = h1[:61]
        h1b = h1[61:]

        # Add current file name to header
        hf  = 'I      BC_GALAXEV --- MODEL PARAMETERS:  Generic file name for model = ' + ff
        hf  = hf + (121 - len(hf))*v + 'I'
        hf  = hf.replace(' ','|')
        hfa = hf[:61]
        hfb = hf[61:]

        # Add today's date to header
        d  = ['Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun']
        m  = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
        dt = datetime.now()
        i  = dt.weekday()
        dt = str(dt)
        hh = '                                           '  # Mon Aug 14 05:28:01 2023'
        hh = hh + d[i] + ' ' + m[int(dt[5:7])-1] + ' '  + dt[8:10] + ' ' + dt[11:19] + ' ' + dt[0:4]
        hh = hh.replace(' ','|')
        hda = hh[:61]
        hdb = hh[61:]

        return [hda, hdb, hfa, hfb, h1a, h1b, h2a, h2b]

    # Build fits compatible header
    io = bt.io
    hh = updhdr()
   #prihdr['header1']  = '|||||||||||||||||||||||||||||||||||||||||||Mon|Aug|14|05:28:0'
   #prihdr['header2']  = '1|2023'
    prihdr['header1']  = hh[0]
    prihdr['header2']  = hh[1]
    prihdr['header3']  = 'I------------------------------------------------------------'
    prihdr['header4']  = '------------------------------------------------------------I'
    prihdr['header5']  = 'I||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||'
    prihdr['header6']  = '||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||I'
   #prihdr['header7']  = 'I||||||BC_GALAXEV|---|MODEL|PARAMETERS:||Generic|file|name|fo'
   #prihdr['header8']  = 'r|model|=|cb2019_zXXX_chab_hr_xmilesi_ssp|||||||||||||||||||I'
    prihdr['header7']  = hh[2]
    prihdr['header8']  = hh[3]
    prihdr['header9']  = 'I||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||'
    prihdr['header10']  = '||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||I'
    prihdr['header11']  = 'I||||||TRACKS:|Padova|2013|+|Marigo|2013,|X=0.704,|Y=0.279,|Z'
    prihdr['header12']  = '=0.017|-|MB|PAGB||||||||||||||||||||||||||||||||||||||||||||I'
    prihdr['header13']  = 'I||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||'
    prihdr['header14']  = '||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||I'
    if io == 5:
        prihdr['header15']  = hh[4]
        prihdr['header16']  = hh[5]
        prihdr['header17']  = hh[6]
        prihdr['header18']  = hh[7]
    else:
       #prihdr['header15']  = 'I||||||S.F.R.:|SSP|=|Zero|Length|Burst|at|t|=|0||||||||||||||'
       #prihdr['header16']  = '||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||I'
        prihdr['header15']  = hh[4]
        prihdr['header16']  = hh[5]
        prihdr['header17']  = 'I||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||'
        prihdr['header18']  = '||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||I'
    prihdr['header19']  = 'I||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||'
    prihdr['header20']  = '||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||I'
    prihdr['header21']  = 'I||||||I.M.F.:|Lognormal|+|power|law|||||||||||||||||||||||||'
    prihdr['header22']  = '||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||I'
    prihdr['header23']  = 'I||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||'
    prihdr['header24']  = '||||||||||mass|in||||||number|||||||||||||||||||||||||||||||I'
    prihdr['header25']  = 'I||||||||||||||||||||||||||||||||||||||||||||||||||from|m||||'
    prihdr['header26']  = '|to|m|||||segment||||of|stars|||||||||||||||||||||||||||||||I'
    prihdr['header27']  = 'I||||||||||||||||||||||||||||||||||||||||lognormal|||0.10||||'
    prihdr['header28']  = '|1.00||||||0.4057||||||1.3076|||||||||||||||||||||||||||||||I'
    prihdr['header29']  = 'I||||||||||||||||||||||||||||||||||||||||||x=1.3|||||1.00|||1'
    prihdr['header30']  = '00.00||||||0.5943||||||0.1827|||||||||||||||||||||||||||||||I'
    prihdr['header31']  = 'I||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||'
    prihdr['header32']  = '||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||I'
    prihdr['header33']  = 'I|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||t'
    prihdr['header34']  = 'otals||||||1.0000||||||1.4903||(|0.6710|Mo/star)||||||||||||I'
    prihdr['header35']  = 'I||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||'
    prihdr['header36']  = '||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||I'
    prihdr['header37']  = 'I||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||'
    prihdr['header38']  = '||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||I'
    prihdr['header39']  = 'I||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||'
    prihdr['header40']  = '||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||I'
    prihdr['header41']  = 'I||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||'
    prihdr['header42']  = '||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||I'
    prihdr['header43']  = 'I||||||LISTED:|Rest|frame|properties|||||||||||||||||||||||||'
    prihdr['header44']  = '||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||I'
    prihdr['header45']  = 'I||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||'
    prihdr['header46']  = '||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||I'
    prihdr['header47']  = 'I||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||'
    prihdr['header48']  = '||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||I'
    prihdr['header49']  = 'I||||||||||||||||||||||||||||||(C)|1995-2023|G.|Bruzual|A.|&|'
    prihdr['header50']  = 'S.|Charlot|-|All|Rights|Reserved||||||||||||||||||||||||||||I'
    prihdr['header51']  = 'I------------------------------------------------------------'
    prihdr['header52']  = '------------------------------------------------------------I'
    prihdr['header53']  = 'END'

    # Write fits file
    f = name
    hdul.writeto(f, overwrite=True)

def getb(a):
    # Check command line argument
    b = ''.join(str(e) for e in a)
    b,k = deca(b)
    return b,k

def getc(a):
    # Ask for time steps to plot
    a = input(a)
    a,k = deca(a)
    return a,k

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
        hd.examples()
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

def nsed(n,a,jfile):
    # Build file name and open output file

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

    # Get requested wavelength range
    w1 = 0.
    w2 = 1.e10
    w0 = 0.
    f0 = 0.
    nu = False
    if a == 0:
        a=''
    if len(a) > 0:
        a = a.replace(',',' ')
        a = a.split()
        w1 = float(a[0])
        if w1 < 0:
            nu = True
            w1 = -w1
        if len(a) > 1:
            w2 = float(a[1])
        if len(a) > 2:
            w0 = float(a[2])
            f0 = 1.
        if len(a) > 3:
            f0 = float(a[3])
    bt.nu = nu

    # get records to list
    i,x = geti(n,t)

    # build file name
    ifile = bt.lfits
    if jfile == '':
        o   = bs.lnam(ifile)
        o   = o.replace('fits',x)
    else:
        o = jfile
    o   = (os.environ.get('glxout') + '/' + o).replace('//','/')
    of  = open(o,'w')
    print('                            Writing records:'  + ''.join('%4i' % i[v] for v in range(1,len(i))))
    print('                              to file: ' + o + '    (will take some time...)')
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

    # Check if Fnu has been requested
    cn = [1.]*(len(w))
    if bt.nu:
        fc= 1.E-8/3.E10
        for l in range(len(w)):
            cn.append(fc*w[l]**2)

    # Check if normalization has been requested
    fn = [1.]*(len(t)+1)
    if w0 > 0:
        jn = np.searchsorted(w, w0)
        for j in range(1,len(i)):
            fn[int(i[j])] = f0/f[jn][int(i[j])] / cn[jn]

    # Write model sed's to file
    for l in range(len(w)):
        if w[l] >= w1 and w[l] <= w2:
            s = f'{w[l]:10.4E}'
            for j in range(1,len(i)):
               #s += f'{f[l][int(i[j])]:11.4E}'
                s += f'{ f[l][int(i[j])] * fn[int(i[j])] * cn[l] :11.4E}'
            of.write(s+'\n')
        elif w[l] > w2:
            break
    of.close()
    print()

def action(n,k):
    # Writes to file required output
    ifile = bt.lfits
    if k == 0:
        hd.ytab(ifile,w,t,e,m,d,p,v,a)
    elif k == -2:
        bcfits2ised(ifile)
    elif k == -3:
        bcfits2ascii(ifile)
    elif k > 0:
       #nsed(n,0,0,'')
        nsed(n,0,'')

def bcf2t(g):
    # Transforms fits file to text tables
    global w,f,t,h,m,d,p,a,s,e,v

    # Execute option selected from the command line

    # Read command line
    k = len(g)

    # Superseded
    #if k < 2:
    #    print('Usage:  bcfits2txt.py model_file.fits')
    #else:
    #    # Read fits file
    #    ifile = g[1]
    if g[0] != bt.lfits:
        bt.lfits = ''
        bt.fw = g[0]

    # Read fits file
    if bt.lfits == '':
        ifile = bt.fw
        ifile = bs.fcheck(ifile)
        w,f,t,e,s,m,p,d,v,a,h = bcfits(ifile)
    else:
        ifile = bt.lfits
        w,f,t,e,s,m,p,d,v,a,h = bt.bc
    if not jupy:
        hd.multhead(1,ifile,t,w,e)
    else:
        hd.multhead(0,ifile,t,w,e)

    # Check for widget command
    if not jupy:
        if k == 1:
            print (' Select records to extract:')
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
                v=' '
                n,k = getc(20*v + 'Choice: ')
                if len(n) <= 0:
                    break
                action(n,k)
        else:
            # command line mode
            n,k = getb(g[1])
            action(n,k)
    else:
        # Check information requested
        if k > 1 and g[1] == 'q':
            hd.multhead(1,ifile,t,w,e)
            return
        if k > 1 and g[1] == 's':
            hd.multhead(2,ifile,t,w,e)
            return
        n,k = deca(str(g[1]))
        action(n,k)

def bcfits2ised(ifile):
    global w,f,t,h,m,d,p,a,s,e,v

    # Writes fortran compatible *.ised2 binary file

    # Open ised2 binary file
    n = bs.lnam(ifile)
    n = n.replace('fits','ised2')
    n = (os.environ.get('glxout') + '/' + n).replace('//','/')
    o = FortranFile(n, 'w')
    print('                            Writing file: ' + n + '    (will take some time...)\n')

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

def bcfits2ascii(ifile):
    global w,f,t,h,m,d,p,a,s,e,v

    # Writes fits file in ASCII format compatible with BC03 release

    # Open ised_ASCII file
    n = ifile.replace('fits','ised_ASCII')
    n = (os.environ.get('glxout') + '/' + n).replace('//','/')
    o = open(n,'w')
    print('                            Writing file: ' + n + '    (will take some time...)\n')

    # Define integer parameters
    nb    = len(t)		# number of time steps
    nw    = len(w)		# number of wavelength points
    io    = s[1][0]		# SFR code
    nskip = s[1][1]		# number of dark SEDs to skip
    kdeff = s[1][2]		# BC/CB model code
    sissa = s[1][3]		#  1 if SISSA 2020 stellar tracks used
    nskip = min(0,nskip)	# Check for cb2022 models which contain dark time steps

    # Write header
    hd.myheader2(o,w,t,e,-1)

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

def glxpl(file1):
    # python version of galaxevpl code
    global w,f,t,e

    # Ask for input parameters to run galaxevpl program
    bs.bchead()
    print(' GALAXEVPL: This program will write a formatted file with galaxy s.e.d.''s')
    print('            read from an unformatted fits table file.')
    print()
    f1 = input(' BC_GALAXEV SSP sed in file [' +     bs.lnam(file1)  + '] = ')
    if len(f1) <= 0:
        f1 = bs.lnam(file1)
    else:
        f1 = bs.zrep(file1,f1)

    # Read BC/CB model in fits file
    w,f,t,e,*nouse = bcfits(f1)
    hd.multhead(1,f1,t,w,e)

    # Ask for time steps to extract
    a = input(' Enter age of up to 250 sed''s = ')
    n, k = deca(a)

    # Ask for wavelength range
    print()
    print(' Enter desired wavelength range = [W1,W2] (default: full range).')
    print('     o If you want all s.e.d.''s scaled to flux = F0 at lambda = W0, enter the')
    print('          desired values (default: no scaling)')
    print('     o If you want the output as Fnu vs. lambda, enter W1 with a minus sign.')
   #print('     o If you want all s.e.d.''s scaled so that F0 is the flux measured through')
   #print('          filter N at redshift Z, enter W0 = -N')
   #a = input('W1,W2,W0,F0,Z = ')
    print()
    a = input(' W1,W2,W0,F0 = ')

    # Ask for output file name
    print()
    b = input(' Output file name [enter for default] = ')
    if len(b) <= 0:
        b = ''
    nsed(n,a,b)

def wglxpl(file1):
    # Runs galaxevpl following data entered in widgets.
    global w,f,t,h,m,d,p,a,s,e,v

    # Read BC/CB model in fits file
    w,f,t,e,*nouse = bcfits(file1)
    hd.multhead(1,file1,t,w,e)

    # Selected records
    rm = bt.rm
    n, k = deca(rm[1])
    nsed(n,rm[2],rm[3])

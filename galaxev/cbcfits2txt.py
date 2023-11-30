#!/usr/bin/env python3

# G. Bruzual (June 2022)

import sys
import math
import numpy as np
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def ctab():
    # Writes header in BC tables
#   a = '#          '
    x = ['.colors', '.mags', '.ABmag', '.lsindx']
    f = lnam(ifile.replace('.fits',''))
    o = '.'
    print('            Creating files: ' + f + '.colors + (.mags, .ABmag, .lsindx)')
    print()
    o = o + '/' + f
    a = m['log_age']				# time scale without pre MS evolution (agrees with t for BC03 and C&B models)
    a[0] = 0.
    for l in range(4):
        of = o + x[l]
        of = open(of,'w')

        if l == 0:
            # Colors
            h = [ 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
            u = [0]*10
            # Header for file .colors
            h1 = '#    (1)         (2)      (3)       (4)       (5)       (6)       (7)       (8)       (9)      (10)      (11)'
            h2 = '#'
            h3 = '#log-age-yr     U-B       B-V       V-R       R-I       I-J       J-H       V-I       V-J       V-H       V-K'
            of.write(h1 + '\n')
            of.write(h2 + '\n')
            of.write(h3 + '\n')
            for r in range(len(a)):
                for j in range(10):
                    u[j] = c[r][h[j]-1]
                s = '%10.6f' % a[r] + ' ' + ' '.join('%9.4f' % v for v in u)
                of.write(s + '\n')
            of.close()

        if l == 1:
            # Vega mags
            h = [ 2, 3, 4, 5, 6, 7, 8, 9]
            u = [0]*8
            # Header for file .mags
            h1 = '#    (1)       (2)       (3)       (4)       (5)       (6)       (7)       (8)       (9)'
            h2 = '#'
            h3 = '#log-age-yr     U         B         V         R         I         J         H         K'
            of.write(h1 + '\n')
            of.write(h2 + '\n')
            of.write(h3 + '\n')
            for r in range(len(a)):
                for j in range(8):
                    u[j] = m[r][h[j]-1]
                s = '%10.6f' % a[r] + ' ' + ' '.join('%9.4f' % v for v in u)
                of.write(s + '\n')
            of.close()

        if l == 2:
            # AB mags
            h = [2, 3, 4, 5, 6]
            u = [0]*5
            # Header for file .ABmag
            h1 = '#    (1)       (2)       (3)       (4)       (5)       (6)'
            h2 = '#'
            h3 = '#log-age-yr    u_AB      g_AB      r_AB      i_AB      z_AB'
            of.write(h1 + '\n')
            of.write(h2 + '\n')
            of.write(h3 + '\n')
            for r in range(len(a)):
                for j in range(5):
                    u[j] = b[r][h[j]-1]
                s = '%10.6f' % a[r] + ' ' + ' '.join('%9.4f' % v for v in u)
                of.write(s + '\n')
            of.close()

        if l == 3:
            # Lick indices
            h = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36]
            u = [0]*35
            # Header for file .lsindx
            h1 = '#     (1)       (2)      (3)      (4)      (5)      (6)      (7)      (8)      (9)     (10)     (11)     (12)     (13)     (14)     (15)     (16)     (17)     (18)     (19)     (20)     (21)     (22)     (23)     (24)     (25)     (26)     (27)     (28)     (29)     (30)     (31)     (32)     (33)     (34)     (35)     (36)'
            h2 = '#'
            h3 = '#log-age-yr     CN1      CN2   Ca4227    G4300   Fe4383   Ca4455   Fe4531    C4668    Hbeta   Fe5015      Mg1      Mg2      Mgb   Fe5270   Fe5335   Fe5406   Fe5709   Fe5782      NaD     TiO1     TiO2  HdeltaA  HgammaA  HdeltaF  HgammaF   Ca8498   Ca8542   Ca8662   Mg8807    B4000  D4000vn  H8-3889  H9-3835 H10-3798     CaHK'
            of.write(h1 + '\n')
            of.write(h2 + '\n')
            of.write(h3 + '\n')
            for r in range(len(a)):
                for j in range(35):
                    u[j] = d[r][h[j]-1]
                s = '%10.6f' % a[r] + ' ' + ' '.join('%8.4f' % v for v in u)
                of.write(s + '\n')
            of.close()

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

def examples():
    # Shows some examples of how to use this procedure
    print ('')
    print ('                  Examples: Screen input                       Command line                                               Output')
    print ("                    Choice: 0.001 0.05 0.1 1 10 11 12 13 14    cbcfits2txt file.fits '0.001 0.05 0.1 1 10 11 12 13 14'    SEDs t = 1Myr, 50 Myr, 100Myr, 1, 10, 11, 12, 13, 14 Gyr")
    print ("                            -1 10 20 50 100 150 175 200 220    cbcfits2txt file.fits '-1 10 20 50 100 150 175 200 220'    SEDs N = 1, 10, 20, 50, 100, 150, 175, 200, 220")
    print ('                            200:210                            cbcfits2txt file.fits 200:210                              SEDs N = 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210')
    print ('                            a                                  cbcfits2txt file.fits a                                    SEDs N = 1 to 221')
    print ('                            t                                  cbcfits2txt file.fits t                                    Tables with model properties')
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

def lnam(s):
    # Returns file name without directory path
    s = s.split('/')
    return s[len(s)-1]

def geta(a):
    # Ask for time steps to plot
    a = input(a)
    a,k = deca(a)
    return a,k

def getb(a):
    # Check command line argument
    b = ''.join(str(e) for e in a)
    b,k = deca(b)
    return b,k

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

def read_cbcfits(ifile,k):
    # Reads fits table with CBC model
    import warnings
    import os
    from   astropy.table import Table
    from   astropy.io import fits
    from   astropy.utils.exceptions import AstropyWarning
    
    # Read file in fits table format
    with fits.open(ifile) as hdul:
        warnings.simplefilter('ignore', category=AstropyWarning)
        # Use Table to read the fits table
        print('Reading file... (will take a few seconds)')
        f = Table.read(hdul,hdu=1)	# hdu = 1 => SED's (luminosity vs. wavelength)
        c = Table.read(hdul,hdu=2)	# hdu = 2 => photometric colors (Johnson's filters)
        m = Table.read(hdul,hdu=3)	# hdu = 3 => photometric Vega magnitudes
        b = Table.read(hdul,hdu=4)	# hdu = 4 => photometric AB magnitudes
        d = Table.read(hdul,hdu=5)	# hdu = 5 => line spectral indices
        w = f['Wavelength']		# wavelength array
        t = 10.**m['log_age']		# time scale in yr
        t[0] = 0.
    if k==0:
        fdoc(w,t)
    return w,f,t,c,m,b,d

def action(n,k):
    # Writes to file required output
    if k == 0:
        ctab()
    elif k == -2:
        bcfits2ised()
    elif (k > 0):
        nsed(n)



# Execute selected option - Command Line version

# Read command line
g = sys.argv
k = len(g)

if k < 2:
    print('Usage:  cbcfits2txt.py model_file.fits')
else:
    # Read fits file
    ifile = g[1]
    w,f,t,c,m,b,d = read_cbcfits(ifile,0)

    if k == 2:
        print ('Select records to extract:')
        print ('                            t1 t2 t3... (age in Gyr of SED''s)')
        print ('                           -N1 N2 N3... (record number of SED''s. Note - sign in front of N1)')
        print ('                            N1:N2       (all SED''s with N1 <= N <= N2)')
        print ('                            a           (all SED''s in fits file)')
        print ('                            t           (write tables with different properties of SSP model)')
        print ('                            h           (enter your choice in the command line, examples)')
        print ('                            ENTER       (exits loop)')
        while True:
            # Ask for record to extract
            n,k = geta('                    Choice: ')
            if len(n) <= 0:
                break
            action(n,k)
    else:
        # command line mode
        n,k = getb(g[2])
        action(n,k)

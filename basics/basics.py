import os
import builtins
import numpy as np
from astropy.table import Table
from termcolor import colored, cprint

def bp(a):
    # Prints text in blue color
    # bp = colored(a, 'blue')
    bp = colored(a, 'green')
    return bp

def gp(a):
    # Prints text in blue color
    gp = colored(a, 'green')
    return gp

def rp(a):
    # Prints text in blue color
    rp = colored(a, 'red')
    return rp

def lfile(z):
    # Read file names from temporary files:
    if z=='-':
        # Use last fits file built by any code
        n = open(os.environ.get('glxtmp') + '/__lastfits__.txt','r')
    elif z=='+':
        # Use indicated file as input for programs run in batch mode
        n = open(os.environ.get('glxtmp') + '/__inputFile__.txt','r')
    elif z=='*':
        # Use indicated file as output for programs run in batch mode
        n = open(os.environ.get('glxtmp') + '/__outputFile__.txt','r')
    f = n.read()
    n.close()
    return f

def fl(f,o):
    if o == 0:
        if len(f) == 1:
            f = [os.environ.get('file00')]
        else:
            f = f[1:]
    else:
        if len(f) == 1:
            f = []
        elif f[1] == '-':
            # Use last fits file built by any code
            f = lfile('-')
        else:
            if f[1].find('z') == 0:
                f[1] = zrep('',str(f[1]))
            f = f[1:]
    return f

def fr(f,o):
    if o == 0:
        if len(f) == 2:
            f = [os.environ.get('file00')]
        else:
            f = f[2:]
    else:
        if len(f) == 2:
            f = []
        elif f[2] == '-':
            # Use last fits file built by any code
            f = lfile('-')
        else:
            if f[2].find('z') == 0:
                f[2] = zrep('',str(f[1]))
            f = f[2:]
    return f

def head():
    print('Galaxy Spectral Evolution Library - pyGALAXEV')
    print('Python Version (C) 2020-2023 - G. Bruzual and S. Charlot - All Rights Reserved')
    print('Last revision: ' + bp(fdate(os.environ.get('glxpyl') + '/.tar.timestamp',1)))

def bchead():
    print()
    print(' Galaxy Spectral Evolution Library (GALAXEV)')
    print(' python version (C) 2019-2023 - G. Bruzual and S. Charlot - All Rights Reserved')
    print()

def fdate(f,i):
    # Returns file creation date
    import datetime
    from   datetime import date
    try:
        mtime = os.path.getmtime(f)
    except OSError:
        mtime = 0
    last_modified_date = datetime.datetime.fromtimestamp(mtime)
    s = str(last_modified_date)
    y = int(s[0:4])
    m = int(s[5:7])
    d = int(s[8:10])
    h = s[11:19]
    o = str(date(y,m,d).ctime())
    if i == 0:
        s = o[0:10] + ' ' + h + ' ' + s[0:4]
    else:
        s = o[0:10] + ' ' + s[0:4]
    return s

def dcode(f,o):
    # Modifies filename according to seleted code 'o' for CB20219 models
    if o == 's':
        f = f.replace('chab','salp')
    elif o == 'k':
        f = f.replace('chab','kroup')
    elif o == 'v':
        f = f.replace('chab','v0p30')
    elif o == '1' and '_v' in f:
        f = f.replace('v0p30','v0p80')
    elif o == '2' and '_v' in f:
        f = f.replace('v0p30','v1p00')
    elif o == '3' and '_v' in f:
        f = f.replace('v0p30','v1p30')
    elif o == '4' and '_v' in f:
        f = f.replace('v0p30','v1p50')
    elif o == '5' and '_v' in f:
        f = f.replace('v0p30','v1p80')
    elif o == '6' and '_v' in f:
        f = f.replace('v0p30','v2p00')
    elif o == '7' and '_v' in f:
        f = f.replace('v0p30','v2p30')
    elif o == '8' and '_v' in f:
        f = f.replace('v0p30','v2p80')
    elif o == '9' and '_v' in f:
        f = f.replace('v0p30','v3p30')
    elif o == '0':
        f = f.replace('_hr','_MU010_hr')
    elif o == '3':
        f = f.replace('_hr','_MU300_hr')
    elif o == '6':
        f = f.replace('_hr','_MU600_hr')
    elif o == 'h':
        f = f.replace('_hr_','_er_')
    elif o == 'i':
        f = f.replace('xmilesi','xindous')
    elif o == 't':
        f = f.replace('xmilesi','stelib')
    elif o == 'b':
        f = f.replace('_hr_xmilesi','_lr_BaSeL')
    return f

def fcheck(file):
    # Checks if fits file exists in various directories
    file = lnam(file).strip()
    if file.find('.fits') < 0:
        file = file + '.fits'
    glxext = sspdir(file)
    while True:
        if len(file) < 5:
            file = ldir(file)
        if os.path.isfile(glxext + file):
            file = glxext + file
           #print('ext','   ',file)
            break
        elif os.path.isfile(os.getcwd() + '/' + file):
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
        else:
            file = input(' File ' + rp(file) + ' does not exist. Enter new file name = ')
    return file

def acheck(file):
    # Checks if any file exists in various directories
    file = lnam(file)
    glxext = sspdir(file)
    while True:
        if os.path.isfile(os.getcwd() + '/' + file):
            file = os.getcwd() + '/' + file
            break
        elif os.path.isfile(os.environ.get('glxssp') + '/' + file):
            file = os.environ.get('glxssp') + '/' + file
            break
        elif os.path.isfile(os.environ.get('glxout') + '/' + file):
            file = os.environ.get('glxout') + '/' + file
            break
        elif os.path.isfile(glxext + file):
            file = glxext + file
            break
        else:
            file = input(' File ' +    file  + ' does not exist. Enter new file name = ')
    return file

def inputc(a):
    # Ask for file name and check that file exists
    file = input(a + ' = ')
    if len(file) > 0:
        if len(file) < 4:
            file = ldir(file)
        file = fcheck(file)
    return file

def lnam(s):
    # Returns file name without directory path
    if isinstance(s, list):
        # s is a list
        s = s[0]
    if '/' in s:     # s is a string
        s = s.split("/")
        return s[len(s)-1]
    return s

def zrep(f,z):
    # Replaces metallicity in file name
    if z == '-' or z== '+' or z== '*':
        f = lfile(z)
        return f
    if z.find('z') != 0:
        return z
    else:
        if f == '':
            f = 'cb2019_z017_chab_hr_xmilesi_ssp.fits'
        i = f.find('_z')
        j = f[i+1:].find('_')
        f = f[:i+1] + z + f[i+j+1:]
        return f

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
#   d = False
    q = True
    m = '/'
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
        m    = amup(filew)
#       d = True
    elif 'bc2022' in filew:
        mdir = 'BC22'
    elif 'cb2022' in filew:
        mdir = 'CB22'
    else:
        q = False
        mdir = ''

    # IMF:
    if q:
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
#   if d:
#       m = amup(filew)
#   else:
#       m = '/'

    if q:
        mdir = os.environ.get('glxssp') + '/' + mdir + m
    else:
        mdir = './'
    return mdir

def ldir(o):
    if (o == 'l'):
        print ()
        print ('To access a CB/BC model you may enter the ' + bp('full file name') +', e.g.:')
        print (bp('   0:') + ' cb2019_z0000_chab_hr_xmilesi_ssp.fits           ' + bp('8:') + ' cb2019_z008_chab_hr_xmilesi_ssp.fits')
        print (bp('   1:') + ' cb2019_z0001_chab_hr_xmilesi_ssp.fits           ' + bp('9:') + ' cb2019_z010_chab_hr_xmilesi_ssp.fits')
        print (bp('   2:') + ' cb2019_z0002_chab_hr_xmilesi_ssp.fits           ' + bp('a:') + ' cb2019_z014_chab_hr_xmilesi_ssp.fits')
        print (bp('   3:') + ' cb2019_z0005_chab_hr_xmilesi_ssp.fits           ' + bp('b:') + ' cb2019_z017_chab_hr_xmilesi_ssp.fits')
        print (bp('   4:') + ' cb2019_z001_chab_hr_xmilesi_ssp.fits            ' + bp('c:') + ' cb2019_z020_chab_hr_xmilesi_ssp.fits')
        print (bp('   5:') + ' cb2019_z002_chab_hr_xmilesi_ssp.fits            ' + bp('d:') + ' cb2019_z030_chab_hr_xmilesi_ssp.fits')
        print (bp('   6:') + ' cb2019_z004_chab_hr_xmilesi_ssp.fits            ' + bp('e:') + ' cb2019_z040_chab_hr_xmilesi_ssp.fits')
        print (bp('   7:') + ' cb2019_z006_chab_hr_xmilesi_ssp.fits            ' + bp('f:') + ' cb2019_z060_chab_hr_xmilesi_ssp.fits')
        print ()
       #print (bp('   a:') + ' cb2019_z014_v0p30_hr_xmilesi_ssp.fits           ' + bp('a:') + ' cb2019_z014_v3p30_hr_xmilesi_ssp.fits')
       #print (bp('   b:') + ' cb2019_z017_v0p30_hr_xmilesi_ssp.fits           ' + bp('b:') + ' cb2019_z017_v3p30_hr_xmilesi_ssp.fits')
       #print (bp('   c:') + ' cb2019_z020_v0p30_hr_xmilesi_ssp.fits           ' + bp('c:') + ' cb2019_z020_v3p30_hr_xmilesi_ssp.fits')
       #print (bp('   d:') + ' cb2019_z030_v0p30_hr_xmilesi_ssp.fits           ' + bp('d:') + ' cb2019_z030_v3p30_hr_xmilesi_ssp.fits')
       #print (bp('   e:') + ' cb2019_z040_v0p30_hr_xmilesi_ssp.fits           ' + bp('e:') + ' cb2019_z040_v3p30_hr_xmilesi_ssp.fits')
       #print ()
        print (bp('   A:') + ' bc2003_z0001_chab_hr_xmilesi_ssp.fits           ' + bp('D:') + ' bc2003_z008_chab_hr_xmilesi_ssp.fits')
        print (bp('   B:') + ' bc2003_z0004_chab_hr_xmilesi_ssp.fits           ' + bp('E:') + ' bc2003_z020_chab_hr_xmilesi_ssp.fits')
        print (bp('   C:') + ' bc2003_z004_chab_hr_xmilesi_ssp.fits            ' + bp('F:') + ' bc2003_z050_chab_hr_xmilesi_ssp.fits')
        print ()
       #print ('or the corresponding ' + bp('zlim') + ' code, defined as follows:')
       #print (bp('    z: ') + '= (0-f)       hexadecimal number indicated above next to the file name. No default value for ' + bp('z') + '.')
       #print (bp('    l: ') + '= (' + bp('m') + ',h,i,t,b) to select the Miles, Miles+, IndoUS, Stelib, or BaSel spectral library')
       #print (bp('    i: ') + '= (' + bp('c') + ',k,s,v)   to select the Chabrier, Kroupa, Salpeter or Vazdekis et al. IMF')
       #print (bp('    m: ') + '= (0,' + bp('1') + ',3,6)   to select Mup = 10, 100, 300, or 600 Mo for the CB2019 Chabrier, Kroupa and Salpeter IMF models')
       #print (bp('       ') + '= (' + bp('1') + ',2,7)     to select the BC2003 (Padova 1994 tracks), CB2003 (Padova 2000 tracks), or CB2007 models')
       #print (bp('       ') + '= (' + bp('0') + ':9)       to select the Vazdekis IMF slope x in the CB2019 models:')
       #print ('                       x = 0.3 (0), 0.8 (1), 1 (2), 1.3 (3), 1.5 (4), 1.8 (5), 2 (6), 2.3 (7), 2.8 (8), 3.3 (9)')
       #print ()
       #print ('                     You must enter your choice for ' + bp('z') + '. The ' + bp('lim ') + 'codes are optional and can be entered in any order.')
       #print ('                     Default values are indicated in ' + bp('blue') + '. Enter ' + bp('?') + ' for examples of ' + bp('zlim') + ' code usage.')
       #print ()
    elif (o=='q'):
        sys.exit()
    elif (o=='h'):
        myhelp()
        print ()
        f = inputc('BC/CB model in fits table file name (model 1)')
        if len(f) <= 0:
            sys.exit()
        return f
    elif (o=='?'):
        print ()
        print (bp('      zlim') + ' code usage:')
        print (bp('         a: ') + 'cb2019_z014_chab_hr_xmilesi_ssp.fits         => Z=0.014, Chabrier IMF, Miles  spectral libray, Mu=100 Mo')
        print (bp('         a: ') + 'cb2019_z014_chab_hr_xmilesi_ssp.fits         => Z=0.014, Chabrier IMF, Miles  spectral libray, Mu=100 Mo')
        print (bp('        as: ') + 'cb2019_z014_salp_hr_xmilesi_ssp.fits         => Z=0.014, Salpeter IMF, Miles  spectral libray, Mu=100 Mo')
        print (bp('        ak: ') + 'cb2019_z014_kroup_hr_xmilesi_ssp.fits        => Z=0.014, Kroupa   IMF, Miles  spectral libray, Mu=100 Mo')
        print (bp('        a3: ') + 'cb2019_z014_chab_MU300_hr_xmilesi_ssp.fits   => Z=0.014, Chabrier IMF, Miles  spectral libray, Mu=300 Mo')
        print (bp('        a6: ') + 'cb2019_z014_chab_MU600_hr_xmilesi_ssp.fits   => Z=0.014, Chabrier IMF, Miles  spectral libray, Mu=600 Mo')
        print (bp('       8bs: ') + 'cb2019_z008_salp_lr_BaSeL_ssp.fits           => Z=0.008, Salpeter IMF, BaSeL  spectral libray, Mu=100 Mo')
        print (bp('      43tk: ') + 'cb2019_z001_kroup_MU300_hr_stelib_ssp.fits   => Z=0.001, Kroupa   IMF, Stelib spectral libray, Mu=300 Mo')
        print (bp('       7v0: ') + 'cb2019_z006_v0p30_hr_xmilesi_ssp.fits        => Z=0.006, Vazdekis IMF, Miles  spectral libray, Mu=100 Mo')
        print (bp('      atv9: ') + 'cb2019_z014_v3p30_hr_stelib_ssp.fits         => Z=0.014, Vazdekis IMF, Stelib spectral libray, Mu=100 Mo')
        print ()
        print ('     Requested files are searched first in the current directory ' + bp('./') + ', then in the distributed SSP model directory ' + bp('$glxssp') +',')
        print ('     and then in the output directory ' + bp('$glxout') + '. Output from the various codes is written to ' + bp('$glxout') + ' and temporary files to ' + bp('$glxtmp') + '.')
        print ()
        # f = inputc('BC/CB model in fits table file name (model 1)')
        # return f
    else:
        # Returns desired model
        c = o[0]
        if c=='A' or c=='B' or c=='C' or c=='D' or c=='E' or c=='F':
            # BC2003 models
            z = [ '', '', '', '', '', '', '', '', '', '', 'z0001', 'z0004', 'z004', 'z008', 'z020', 'z050' ]
            f = 'bc2003_' + z[int(o[0],16)] + '_chab_hr_xmilesi_ssp.fits'
            if len(o) > 1:
                f = ecode(f,o[1])
            if len(o) > 2:
                f = ecode(f,o[2])
            if len(o) > 3:
                f = ecode(f,o[3])
            return f
        elif c=='0' or c=='1' or c=='2' or c=='3' or c=='4' or c=='5' or c=='6' or c=='7' or c=='8' or c=='9' or c=='a' or c=='b' or c=='c' or c=='d' or c=='e' or c=='f':
            # CB2019 model
            z = [ 'z0000', 'z0001', 'z0002', 'z0005', 'z001', 'z002', 'z004', 'z006', 'z008', 'z010', 'z014', 'z017', 'z020', 'z030', 'z040', 'z060' ]
            f = 'cb2019_' + z[int(o[0],16)] + '_chab_hr_xmilesi_ssp.fits'
            if len(o) > 1:
                f = dcode(f,o[1])
            if len(o) > 2:
                f = dcode(f,o[2])
            if len(o) > 3:
                f = dcode(f,o[3])
            return f
        else:
            return o

def trapz2(x,y,x1,x2):
    # Finds area below y(x) from x=X1 to x =X2, X1 and X2 are not necessarily values of x(i) but must fulfill X1>=x(1) and X2<=x(n)

    # Check for values outside valid range
    #if x1 < x[0]:
    #    print('TRAPZ2: X1 outside valid range:',x1,' using x1 =',x[0])
    #    x1 = x[0]
    #if x2 > x[len(x)]:
    #    print('TRAPZ2: X2 outside valid range:',x2,' using x2 =',x[len(x)-1])
    #    x2 = x[len(x)]-1
    #    quit()
    x1 = max(x1,x[0])
    x2 = min(x2,x[len(x)-1])

    # Locate x1, x2 in array x
    i1, i2 = np.searchsorted(x,[x1,x2])
    # Interpolate (x,y) at (x1,x2)
    y1 = np.interp(x1,x,y)
    y2 = np.interp(x2,x,y)
    # Build integration array
    u=[]
    v=[]
    u.append(x1)
    v.append(y1)
    for i in range(i1,i2):
        u.append(x[i])
        v.append(y[i])
    u.append(x2)
    v.append(y2)
    # Compute area below curve with trapezoidal rule
    a = np.trapz(v,u)
    return a

def myhelp():
    # Prints a brief summary of available tasks
    s = '         '
    print ()
    print ('   The following tasks are available to explore the GALAXEV models:')
    print ()
    print (s + bp('ascii_tables') + ' - writes tables with the properties of the selected model displayed with options ' + bp('s,n,m,c,i,p'))
    print (s + bp('model_ID') + '     - displays the basic ingredients defining the selected model')
    print (s + bp('csp_galaxev') + '  - creates a composite stellar population model for the indicated SFR')
    print (s + bp('add_bursts ') + '  - adds two stellar population models separated in time and in proportions indicated by the user')
    print (s + bp('cm_evolution')+  ' - computes the observer frame apparent magnitude and colors versus redshift in the selected filter(s)')
    print (s + bp('zmag        ')+  ' - returns the observer frame apparent magnitude of the model stellar population at the indicated redshift')
    print (s + bp('galaxevpl   ')+  ' - writes table with the model spectra at selected ages (also possible with options ' + bp('s,n') + ')')
    print (s + bp('fits2ised   ')+  ' - writes part of the content of the fits table file in the ' + bp('.ised') + ' format introduced by ' + bp('Bruzual & Charlot (2003)'))
    print ()
    print (s + bp('File names:')+  '  - Starting with the present release, each GALAXEV SSP model is provided as a single fits table file, for example,')
    print ()
    print (2*s + '                        ' + bp('cb2019_z014_chab_hr_xmilesi_ssp.fits'))
    print (2*s + '                        ' + bp('cb2019_z017_kroup_MU300_er_xmilesi_ssp.fits'))
    print (2*s + '                        ' + bp('cb2019_z004_salp_MU600_lr_BaSeL_ssp.fits'))
    print (2*s + '                        ' + bp('cb2007_z008_chab_hr_xindous_ssp.fits'))
    print (2*s + '                        ' + bp('bc2003_z020_chab_hr_stelib_ssp.fits'))
    print ()
    print ('   where:')
    print (s + bp('cb2019, cb2007, bc2003') + ' indicate the model release,')
    print (s + bp('z014, z017, z004, z008, z020') + ' is the model metallicity ' + bp('Z')+' = 0.014, 0.017, 0.004, 0.008, 0.02 in the examples,')
    print (s + bp('chab, kroup, salp') + ' denote the model IMF: Chabrier, Kroupa, Salpeter,')
    print (s + bp('miles, indous, stelib, BaSeL') + ' indicate the stellar library used in the visible range, and ' + bp('er, hr, lr ') + 'its resolution')
    print (s + bp('MU100, MU300, MU600') + ' indicate the IMF upper mass limit Mu = 100, 300, 600 Mo. If not included, Mu = 100 Mo.')
    print (s + bp('ssp') + ' specifies that the file content corresponds to a simple stellar population (SSP).')
    print (s + 'For details see the ' + bp('pyGALAXEV Release Notes') + '(in preparation).')
    #rint ()
    #rint ()
    #rint ()

def flog(glxdir,name):
    # Shows in screen content of file "name"
    name = os.environ.get(glxdir) + '/' + name
    if 'filters.list' in name:
        print ('List of filters stored in file: ' + os.environ.get('FILTERF').replace('//','/'))
        print()
    with open(name, 'r') as f:
        line = f.read()
        print(line)


import math
import numpy    as np

# Concordance cosmological model.
# The term ‘concordance model’ is used in cosmology to indicate the currently accepted
# and most commonly used cosmological model. It is important to identify a concordance
# model because the measurement of many astrophysical quantities (e.g. distance, radius,
# luminosity and surface brightness) depend upon the cosmological model used. Consequently,
# for ease of comparison if nothing else, the models assumed in different studies should
# at least be similar, if not identical.

# Currently, the concordance model is the Lambda CDM model (which includes cold dark matter
# and a cosmological constant). In this model the Universe is 13.7 billion years old and made
# up of 4% baryonic matter, 23% dark matter and 73% dark energy. The Hubble constant for this
# model is 71 km/s/Mpc and the density of the Universe is very close to the critical value for
# re-collapse. These values were derived from WMAP (Wilkinson Microwave Anisotropy Probe)
# satellite observations of the cosmic microwave background radiation.

# Parameters:

# h = 71
# omega = 0.27
# omega_lambda = 0.73
# tu = 13.7 Gyr

def cosmol():
    global h, omega, omega_lambda, clambda, q, tu, ttg, zf
    # Define default value of cosmological parameters
    h = 71.  ; omega = 0.27 ; omega_lambda = 0.73 ; ttg = 13.5
    clambda,q = cosmol_c(h,omega,omega_lambda)
    tu = tuniverse(h,q,0.,clambda)
    zf = zx(ttg,h,q,clambda)
    # Return cosmological parameters
    return h, q, clambda, tu, ttg, zf, omega, omega_lambda

def cosmol_c(h,omega,omega_lambda):
    # Returns cosmological constant = cosmol_c and parameter q

    #  omega is entered by the user
    #  omega=1.-omega_lambda

    #  cosmological constant
    cosmol_c=omega_lambda/(3.*h**2)

    # compute q = q0 (deceleration parameter)
    if omega_lambda == 0.:
        q = omega/2.
    else:
        q = (3.*omega/2.) - 1.
    return cosmol_c, q

def dismod(h,q,z):
    # Returns cosmological distance modulus

    # h  = Ho in km/sec/Mpc
    # q  = qo
    # dl = Luminosity distance in Mpc

    dismod = 5*np.log10(dl(h,q,z)*1.E6/10.)
    return dismod

def dl(h,q,z):
    # Computes luminosity distance corresponding to a redshift z.
    # Uses Mattig formulae for qo both 0 and non 0
    # Revised January 1991 to implement cosmological constant
    # h = Ho in km/sec/Mpc
    # ******DL is in Mpc******
    global omega0

    def funl(x):
        # For non-zero cosmological constant
        omegainv = 1. / omega0
        funl = 1. / math.sqrt(((x ** 3.) + omegainv) - 1.)
        return funl

    if z <= 0.:
        # 10 pc
        dl=1.e-5
        return dl
    if q == 0:
        dl = ((3.e5 * z) * (1 + (z / 2.))) / h
    elif q > 0.:
        d1 = (q * z) + ((q - 1.) * (math.sqrt(1. + ((2. * q) * z)) - 1.))
        d2 = ((h * q) * q) / 3.e5
        dl = d1 / d2
    elif q < 0.:
        omega0 = (2. * (q + 1.)) / 3.
        aa = 1.
        bb = 1. + z
        s0 = 1.e-10
        s  = 0.
        npts=0
        while True:
            npts=npts+1
            s =  midpnt(funl,aa,bb,s,npts)
            epsr=abs(s-s0)/s0
            if epsr < 1.e-4:
                break
            else:
                s0=s
        dd1 = s
        dd2 = (3.e5 * (1. + z)) / (h * math.sqrt(omega0))
        dl  = dd1 * dd2
    return dl

def ltt(h,q,zi,lamb,age):
    # Define light travel time
    ltt = age - tuniverse(h,q,zi,lamb)
    return ltt

def midpnt(f,a,b,s,n):
    if n == 1:
        s=(b-a)*f(0.5*(a+b))
    else:
        it=3**(n-2)
        tnm=it
        delt=(b-a)/(3.*tnm)
        ddel=delt+delt
        x=a+0.5*delt
        tot=0.
        for j in range(1,it+1):		 #do 11 j=1,it
            tot=tot+f(x)
            x=x+ddel
            tot=tot+f(x)
            x=x+delt
        s=(s+(b-a)*tot/tnm)/3.
    return s

def tuniverse(h, q, z, lamb):
    # Returns age of universe at redshift z
    #  H = Ho in km/sec/Mpc
    #  Q = qo  (if problems with qo = 0, try 0.0001)
    global omega0

    def a(q,z):
        a = (math.sqrt(1. + ((2. * q) * z)) / (1. - (2. * q))) / (1. + z)
        return a

    def b(q):
        b = q / (abs((2. * q) - 1.) ** 1.5)
        return b
    
    def c(q,z):
        c = ((1. - (q * (1. - z))) / q) / (1. + z)
        return c

    def funq(x):
        # For non-zero cosmological constant
        omegainv = 1. / omega0
        funq = math.sqrt(omegainv) / (x*math.sqrt((omegainv-1.)+(1./(x**3.))))
        return funq

    def cosh(x):
        #cosh = alog(x + math.sqrt((x ** 2) - 1.))
        cosh = math.log(x + math.sqrt((x ** 2) - 1.))
        return cosh

    hh0 = h * 0.001022		# h in (billion years)**(-1)

    if lamb != 0.0:
        omega0 = (2. * (q + 1.)) / 3.
        s  = 0.
        aa = 0.
        bb = 1. / (1. + z)
        s0 = 1.e-10
        npts = 0
        while True:
            npts=npts+1
            s = midpnt(funq,aa,bb,s,npts)
            epsr=abs(s-s0)/s0
            if epsr < 1.e-4:
                break
            else:
                s0=s
        t = s
    elif q == 0.0:
        t = 1. / (1. + z)
    elif q < 0.5:
        t = a(q,z) - (b(q) * cosh(c(q,z)))
    elif q == 0.5:
        t = (2. / 3.) / ((1. + z) ** 1.5)
    else:
        t = a(q,z) + (b(q) * cos(c(q,z)))
    tuniverse = t / hh0
    return tuniverse
def zx(tx, h, q, lamb):

    # Returns the value of Z = ZX (redshift) corresponding to a given
    # light travel time TX (measured in Gyr).

    #   H = Ho in km/sec/Mpc
    #   Q = qo (if problems with qo = 0, try 0.0001)

    # Define list of z values
    # data z /0.,.001,.002,.004,.006,.008,.01,.02,.04,.06,.08,.1,.2,.3,.4,.5,.545,.6,.7,.8,.9,.945,1.,1.2,1.4,1.6,1.8,2.,3.,4.,5.,6.,7.,8.,9.,10.,12.,14.,16.,18.,20.,30.,40.,60.,80.,100.,1000./
    zl = ["0.", ".001", ".002", ".004", ".006", ".008", ".01", ".02", ".04", ".06", ".08", ".1", ".2", ".3", ".4", ".5", ".545",
          ".6", ".7", ".8", ".9", ".945", "1.", "1.2", "1.4", "1.6", "1.8", "2.", "3.", "4.", "5.", "6.", "7.", "8.", "9.", "10.",
          "12.", "14.", "16.", "18.", "20.", "30.", "40.", "60.", "80.", "100.", "1000."]
    z = []
    for item in zl:
        z.append(float(item))

    # check for zero tx value
    if tx == 0.:
        zx = 0.
        return zx

    # compute omega
    if lamb != 0.0:
        omega0 = (2. * (q + 1.)) / 3.
        omegainv = 1. / omega0

    # check for maximum age of universe
    age = tuniverse(h,q,0.,lamb)
    if tx >= age:
        zx = -2.
        return zx

    # express Ho in billion years ** (-1)
    hh0 = h * 0.001022

    # check for q=0.5 case
    if q == 0.5:
        zx = ((1. - (((3. * hh0) * tx) / 2.)) ** (- (2. / 3.))) - 1.
        return zx

    # general case
    for j in range (len(z)):		# do j = 1, 47
        if tx <= ltt(h,q,z[j],lamb,age):
            if j > 1:
                zx = z[j-1]
            else:
                zx = z[1]
            zi=zx

            # Iterate to find best value
            for j in range(100000):	# do j = 1, 100000
                if q >= 0.:
                    # zero cosmological constant
                    dz = ((hh0 * ((tx - tuniverse(h,q,0.,lamb)) + tuniverse(h,q,zx,lamb))) * ((1. + zx) ** 2)) * ((1. + ((2. * q) * zx)) ** 0.5)
                else:
                    # non-zero cosmological constant
                    dz = ((hh0 * ((tx - tuniverse(h,q,0.,lamb)) + tuniverse(h,q,zx,lamb))) * (1. + zx)) * math.sqrt(omega0)
                    dz = dz * math.sqrt((((1. + zx) ** 3.) + omegainv) - 1.)
                if abs(dz/zx) < 0.00003:
                    return zx
                else:
                    zx = zx + dz
            return zx


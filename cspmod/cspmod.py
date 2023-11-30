import os
import math
import builtins      as bt
import numpy         as np
import bcfits.bcfits as bc
import basics.basics as bs
import bcfilt.bcfilt as fl
import common.common as cn
import cosmol.cosmol as cl

def dustpar():
    ad = input(' Include attenuation by dust? Y/[N]: ')
    v=' '
    if ad == 'y' or ad == 'Y':
       #print()
        print(13*v + 'Using simple 2-component model of Charlot & Fall (2000)]')
        tv = input(13*v + 'Enter total effective attenuation optical depth: tau_V [1.0] = ')
        mu = input(13*v + 'Enter fraction of tau_V arising from the ambient ISM: mu [0.3] = ')
        if len(tv) <= 0:
            tv='1.0'
        if len(mu) <= 0:
            mu='0.3'
        ad = True
    else:
        tv = '0'
        mu = '0'
        ad = False
    return ad,tv,mu

def tausfr():
    # Exponential SFR (enter tau)
    v=' '
    tau = input(13*v + 'Exponential with e-folding time TAU (Gyr) = ')
    ans = input(13*v + 'Recycle gas ejected by stars Y/[N] = ')
    if tau=='':
        tau = 1.
    else:
        tau = float(tau)
    tmu = 1.-math.exp(-1./tau)
    tau = 1.e9*tau
    if ans == 'Y' or ans == 'y':
        io = 1
        epsilon = get_epsilon()
        ans = 'y'
    else:
        ans = 'n'
        epsilon = '0'
        io = 4
    hdr = [tmu, tau, epsilon]
    return io,tau,ans,epsilon,hdr

def musfr():
    # Exponential SFR (enter mu)
    v=' '
    tmu = input(13*v + 'Exponential with mu_SFR parameter = ')
    ans = input(13*v + 'Recycle gas ejected by stars Y/[N] = ')
    if ans == '':
        ans = 'n'
    if tmu == '':
        tmu = 0.5
    else:
        tmu = float(tmu)
    tau = -1.e9/np.log(1.-tmu)
    if ans == 'Y' or ans == 'y':
        io = 1
        epsilon = get_epsilon()
        ans = 'y'
    else:
        epsilon = '0'
        io = 4
    hdr = [tmu, tau, epsilon]
    return io,tau,ans,epsilon,hdr

def bsfr():
    # Long Burst SFR (enter duration of burst)
    v=' '
    tb = input(13*v + 'Duration of burst (Gyr) = ')
    if tb=='':
        tau = 1.e9
    else:
        tau=1.e9*float(tb)
    tcut=tau
    hdr = [tau]
    return tau,tcut,hdr

def csfr():
    # Constant SFR (enter constant value)
    v=' '
    c = input(13*v + 'Enter SFR in Mo/yr [1] = ')
    if c=='':
        c=1.
    else:
        c = float(c)
    hdr = [c]
    return c,hdr

def dsfr():
    # Delayed SFR (enter constant value)
    v=' '
    print(13*v + 'Delayed SFR as defined in Bruzual (1983)')
    td = input(13*v + 'Maximum in SFR at time TAU (Gyr) = ')
    if td=='':
        td=1.E9
    else:
        td=1.e9*float(td)
    hdr = [td]
    return td,hdr

def lsfr():
    # Linear SFR (enter constant value)
    v=' '
    print (13*v + 'Linearly decreasing SFR')
    tl = input(13*v + 'SFR = 0 at time TAU (Gyr) = ')
    if tl=='':
        tl = 1.E9
    else:
        tl=1.e9*float(tl)
    hdr = [tl]
    return tl,hdr

def usrsfr(t,a):
    # SFR(t) read from 2 column from ASCII file
    global time, usfr, tcut, hdr
    v=' '
    if t < 0:
        if a == '':
            print()
            print(13*v + 'SFR(t) read from disk file.')
            print(13*v + 't(yr), SFR(t) in Mo/yr read from 2 column ascii file.')
            print(13*v + 'Linear interpolation between listed points.')
           #print(13*v + 'Skip lines starting with #')
            a = input(13*v + 'SFR in file name = ')
            bt.hdr = [a]
        a = bs.acheck(a)
        d = bc.nfile(a,1)
        time = []
        usfr = []
        for i in range(len(d)):
            time.append(d[i][0])
            usfr.append(d[i][1])
       #print(' ',len(time),'data points read from file',a)
        # If needed, add a first point at age = 0, and a second point just before first usr point.
        if time[0] > 0:
            t1 = 0
            t2 = 0.99995*time[0]
            u1 = 0
            u2 = 0
            time.insert(0,t2)
            time.insert(0,t1)
            usfr.insert(0,u2)
            usfr.insert(0,u1)
        tcut = time[len(time)-1]
        print()
        print(15*v,len(time),'data points in user defined SFR in file',a)

    # Interpolate SFR at time t
    u = np.interp(t, time, usfr, left=0, right=0)
    return u

def chensfr():
    # Double exponential SFR after Chen et al. (2012)
    global ampr, tauj, tau1, tau2, con1, con2

    # Note on the parameters for the Chen et al. SFR as used in the SSAG (Jan-23-2015):
    #   Ayer Gladis y yo encontramos una inconsistencia entre la SFR que calculas con
    #   el programa csp_galaxev para correr las SFH a un z dado y los parámetros de la
    #   SSAG. El problema parece ser la interpretación del tiempo de truncado (tcut,
    #   cuando comienza a decaer más rápidamente la SFR). Ivan calcula estas cosas en
    #   look-back-time, lo que quiere decir que si él dice que una galaxia que se formó
    #   hace 8 Gaños tiene un tiempo de truncado de 1 Gaño, quiere decir que comenzó a
    #   truncarse hace 1 Gaño a partir del presente. Mientras en tus SFR aparecería el
    #   tiempo de truncado como 8-1 G años. Afortunadamente esto solo afecta a las
    #   galaxias con truncado, solo 6 (+1 a z>0) en la muestra y puedo continuar sin
    #   estas galaxias. Alfredo and Gladis.
    print()

    # Express look-back-time Tform as age of universe when galaxy forms
    h = 71. ; omega = 0.27 ; omega_lambda = 0.73
    # Obtain cosmological constant and parameter q
    clambda,q = cl.cosmol_c(h,omega,omega_lambda)
    # age of universe at z = 0
    tuni = cl.tuniverse(h,q,0.,clambda)*1.E9

    def tinp(a):
        # Input value of give time and transforms to yr
        t = input(a)
        if len(t) <= 0:
            t = 0.
        else:
            t = float(t)*1.E9
        return t

    # Ask for parameters for two exponential SFR's
    v=' '
    tauf = tinp(22*v + 'Assume that galaxy formed at look-back-time Tform (Gyr) = ')        # = galaxy age
    tau1 = tinp(17*v + 'For t <= Tj, exponential SFR with e-folding time TAU_1 (Gyr) = ')
    tau2 = tinp(17*v + 'For t >  Tj, exponential SFR with e-folding time TAU_2 (Gyr) = ')
    con1 = 1./tau1
    if tau2 != 0:
        tauj = tinp(22*v + "Join continuously both SFR's at look-back-time Tj (Gyr) = ")
        # Express look-back-time Tj as universe age (from t = 0)
        zauj = cl.zx(tauj*1.E-9,h,q,clambda)
        tauj = tauf - tauj
        # Normalization and continuity of exponential IMF's
        if tauj/tau2 < 100:
            con2 = con1 * math.exp(-tauj/tau1) / math.exp(-tauj/tau2)
        else:
            con2 = 0
    else:
        tau2 = tau1
        con2 = con1
        tauj = 0
    tfrm = tuni - tauf
    zfrm = cl.zx(tfrm*1.E-9,h,q,clambda)

    # Ask for parameters for optional burst of star formation
    ampl = input(11*v + 'Add burst of fractional amplitude A = mass in burst/subyacent mass = ')
    if len(ampl) <= 0:
        ampl = 0
        taui = 0.
        tdur = 0.
    else:
        ampl = float(ampl)
        taui = tinp(34*v + 'Burst starts at look-back-time Tburst (Gyr) = ')
        tdur = tinp(54*v + 'Duration of burst (Gyr) = ')
        # Express look-back-time Tburst as universe age (from t = 0)
        taui = tauf - taui
        zaui = cl.zx(taui*1.E-9,h,q,clambda)
        tend = taui + tdur
        zend = cl.zx((tauf-tend)*1.E-9,h,q,clambda)
    ampr=ampl

    # Build time scale to use in SFR
    time = []
    for i in range(len(t)):
        time.append(t[i])
    # Add steps around tauj if needed
    if tauj > 0:
        time.append(0.9999*tauj)
        time.append(0.99999*tauj)
        time.append(tauj)
        time.append(1.0001*tauj)
        time.append(1.00001*tauj)
    # Add steps around taui if needed
    if ampr > 0:
        time.append(0.9999*taui)
        time.append(0.99999*taui)
        time.append(taui)
        time.append(1.0001*taui)
        time.append(1.00001*taui)
        time.append(0.9999*(taui+tdur))
        time.append(0.99999*(taui+tdur))
        time.append(taui+tdur)
        time.append(1.0001*(taui+tdur))
        time.append(1.00001*(taui+tdur))
        nbr=25
        dtbr=tdur/nbr
        for i in range(nbr+2):
            time.append(taui+i*dtbr)
    # Sort and suppress duplicate entries in array time
    time = np.unique(time)  # unique sorts array time and suppresses duplicate entries (no need for time = np.sort(np.unique(time)))

    # Compute SFR for both exponentials plus burst. Fill in SFR array
    usfr = []
    for i in range(len(time)):
        if time[i] <= tauj:
            u = con1*math.exp(-time[i]/tau1)
        else:
            u = con2*math.exp(-time[i]/tau2)
        usfr.append(u)

    # Compute mass of subyacent population
    tsubya = np.trapz(usfr,time)
    # Scale burst amplitude to obtain required fraction
    # Compute mass of stars formed in burst
    if ampr > 0:
        ampl   = ampr*tsubya/tdur
        tburst = ampl*tdur
        # Add burst population
        for i in range(len(time)):
            if time[i] >= taui and time[i] <= taui+tdur:
                usfr[i] = usfr[i] + ampl

    # Compute total mass in stars at t = tform
    tstelm = bs.trapz2(time,usfr,0.,tauf)
    # Compute mass of stars formed during last Gyr
    tlast1 = bs.trapz2(time,usfr,tauf-1.e9,tauf)

    o = open(os.environ.get('glxout') + '/' + 'py.500', 'w')
    o.write('# SFR parameters after Chen et al. (2012, MNRAS, 421, 314)\n')
    v = ' '
    o.write('#\n')
    s = '# Tuniv  =' + f'{tuni*1.E-9:11.5f} Gyr,          at z = {0.00:6.3f} <-- look-back-time = {0.00*1.E-9:5.2f} Gyr (= age of the universe)' ; o.write(s+'\n')
    s = '# Tform  =' + f'{tauf*1.E-9:11.5f} Gyr,          at z = {zfrm:6.3f} <-- look-back-time = {tauf*1.E-9:5.2f} Gyr (= age of the galaxy)'  ; o.write(s+'\n')
    s = '#' + 31*v + '--> Add ' + f'{tfrm*1.e-9:6.3f} Gyr to model galaxy time scale so that galaxy formation starts {tauf*1.e-9:5.2f} Gyr ago at z = {zfrm:6.3f}' ; o.write(s+'\n')
    o.write('#\n') ; # print(14*v,s )
    if tauj > 0:
        s = '# Tj     =' + f'{tauj*1.E-9:11.5f} Gyr,          at z = {zauj:6.3f} <-- look-back-time = {(tauf-tauj)*1.E-9:5.2f} Gyr' ; o.write(s+'\n')
    if tau1 < 1.E11:
        s ='# TAU_1  =' + f'{tau1*1.E-9:11.5f}' + ' Gyr'       ; o.write(s+'\n')
    else:
        s ='# TAU_1  =' + f'{tau1*1.E-9:10.1f}' + ' Gyr'       ; o.write(s+'\n')
    if tau2 < 1.E11:
        s ='# TAU_2  =' + f'{tau2*1.E-9:11.5f}' + ' Gyr'       ; o.write(s+'\n')
    else:
        s ='# TAU_2  =' + f'{tau2*1.E-9:10.1f}' + ' Gyr'       ; o.write(s+'\n')
    o.write('#\n')
    s = '# Burst amplitude =' + f'{ampr:11.5f}'                ; o.write(s+'\n')
    if ampr > 0:
        s = '# Burst starts at =' + f'{taui*1.E-9:11.5f} Gyr, at z = {zaui:6.3f} <-- look-back-time = {(tauf-taui)*1.E-9:5.2f} Gyr' ; o.write(s+'\n')
        s = '# Burst ends at   =' + f'{tend*1.E-9:11.5f} Gyr, at z = {zend:6.3f} <-- look-back-time = {(tauf-tend)*1.E-9:5.2f} Gyr' ; o.write(s+'\n')
    o.write('#\n')

    # Report stellar mass
    s = 'Mass computed from this SFR:' ; o.write('# ' + s + '\n') # ; print(14*v,s)
    o.write('#\n') ; # print()
    s = 'Mass of stars in subyacent population =' + f'{tsubya:13.5E} Mo' ; o.write('# ' + s + '\n') # ; print(14*v,s)
    if ampr > 0:
        s = 'Mass of stars formed in burst         =' + f'{tburst:13.5E} Mo' ; o.write('# ' + s + '\n') # ; print(14*v,s)
    else:
        tstelm = tsubya
    s = 'Mass of stars in total population     =' + f'{tstelm:13.5E} Mo'     ; o.write('# ' + s + '\n') # ; print(14*v,s)
    s = 'Mass of stars formed during last Gyr  =' + f'{tlast1:13.5E} Mo'     ; o.write('# ' + s + '\n') # ; print(14*v,s)
    if ampr > 0:
        s0 = 'Mass ratio burst/subyacent population =' + f'{ampr:9.5f} = requested burst amplitude'
        s1 = 'Mass ratio burst/subyacent population =' + f'{tburst/(tstelm-tburst):9.5f} = resulting burst amplituder, difference with requested amplitude =' \
            f'{(tburst/(tstelm-tburst) -ampr)/ampr*100.:6.3f}%'
        o.write('# ' + s0 + '\n') # ; print(14*v,s0)
        o.write('# ' + s1 + '\n') # ; print(14*v,s1)
    o.write('#\n') ; # print()

    # Write SFR to table
    o.write('#    t(yr)           SFR (Mo/yr)\n')
    for i in range(len(time)):
        if time[i] <= tauf:
            o.write(str(time[i]) + ' ' + str(usfr[i]) + '\n')
    o.close()

def pgas(tx):
    # Returns amount of processed gas at time tx
    # Interpolate array (so,to) if possible
    lgas = len(to)-1
    if lgas < 0:
        pgas = 0
    elif tx < to[lgas]:
        pgas = np.interp(tx, to, so, left=0, right=0)
    else:
        pgas = so[lgas]
    return pgas

def get_epsilon():
    # Gets fraction of ejected gas to be recycled in stars = epsilon
    v=' '
    print (15*v + 'Epsilon = fraction of ejected gas to be recycled in stars')
    print (15*v + 'Values from 0.001 to 1 have been explored')
    print (15*v + 'Epsilon = 0.001 reproduces old_galaxev mu_SFR=0.50 model')
    print (15*v + 'Epsilon > 1 emulates gas infall')
    print (15*v + 'Epsilon < 0 emulates galactic wind')
    epsilon = input (15*v + 'Epsilon = ')
    if epsilon == '':
        epsilon = 0.
    else:
        epsilon = float(epsilon)
    return epsilon

def cspwidgtparam(file1):
    # Rebuilds csp_galaxev input parameters from data gathered by widget
    global io, tau, ans, epsilon, tcut, tmax, ad, tv, mu, zu, hdr
    tmax  = 20.E9
    i     = file1
    f1    = i[0]
    io    = int(i[1])
    if io == 7:
        u = usrsfr(-1.,i[2])
    else:
        tau = float(i[2])
        if io == 1:
            tmu = 1.-math.exp(-1./tau)
            io = 4 # taumodels with no recycling form widget
        if io !=3:
            tau = tau*1.e9
    tcut  = float(i[3])
    if tcut > 0:
        tcut = tcut*1.E9
    else:
        tcut = tmax
    ans   = i[4]
    ad    = i[5]
    tv    = float(i[6])
    mu    = float(i[7])
    zu    = float(i[8])
    fo    = i[9]
    epsilon = '0'
    bt.io = io
    if io == 1:
        bt.hdr = [tmu, tau, epsilon]
    elif io == 2:
        tcut=tau
        bt.hdr = [tau]
    elif io==3:
        bt.hdr = [tau]
    elif io==4:
        bt.hdr = [tmu, tau, epsilon]
    elif io == 6:
        td=tau
        bt.hdr = [td]
    elif io == 8:
        tl=tau
        bt.hdr = [tl]
    elif io == 7:
        bt.hdr = [i[2]]
    bt.zu = zu
    return f1, fo

def cspparams():
    # SFR parameters for csp_galaxev
    global io, tau, ans, epsilon, tcut, tmax, ad, tv, mu, zu, hdr

    def taustr(t):
        # Return t in Gyr as a string
        t = t*1.E-9
        a = math.modf(t)
        if a[0] == 0:
            s = str(int(a[1]))
        else:
            s = f'{t:6.3f}'.replace(' ','')
        return s

    # Select SFR
    io = sfrpar()
    if io == 7:
        io=9
    elif io == 6:
        io=7
    elif io == 5:
        io=8
    elif io == 4:
        io=6

    # Select SFR parameters for option io
    tmax = 20.E9
    tcut = tmax
    v=' '
    if io == 1:
        # Exponential SFR (enter tau)
        io,tau,ans,epsilon,bt.hdr = tausfr()
        s = 'tau' + taustr(tau) + 'Gyr'
    elif io == -1:
        # Exponential SFR (enter mu)
        io,tau,ans,epsilon,bt.hdr = musfr()
        s = 'tau' + taustr(tau) + 'Gyr'
    elif io == 2:
        # Long Burst SFR (c-model, enter duration of burst)
        tau,tcut,bt.hdr = bsfr()
        s = 'cmod' + taustr(tau) + 'Gyr'
    elif io == 3:
        # Constant SFR (enter constant value)
        tau,bt.hdr = csfr()
        s = 'cons' + taustr(tau*1.E9) + 'Moyr'
    elif io == 6:
        # Delayed SFR (enter time at which maximum SFR occurs)
        tau,bt.hdr = dsfr()
        s = 'dlyd' + taustr(tau) + 'Gyr'
    elif io == 7:
        # Function SFR(t) read 2 column from ASCII file
        u = usrsfr(-1,'')
        s = 'numr'
    elif io == 8:
        # Linear SFR (enter time at which SFR  = 0)
        tau,bt.hdr = lsfr()
        s = 'linr' + taustr(tau) + 'Gyr'
    elif io == 9:
        # Chen et al. double exponential + burst SFR
        chensfr()
        s = ''
        if ampr > 0:
            u = usrsfr(-1,'py.500')
    if io != 2 and io != 7 and io != 9:
        tc = input(13*v + 'Quenching: make SFR = 0 at time > TCUT [' + str(int(tmax*1.E-9)) + ' Gyr] = ')
        if tc != '':
            tcut = 1.e9*float(tc)
    bt.io = io
    bt.sn = s
    print()

    # Dust parameters
    ad,tv,mu = dustpar()
    print()

    # Flux weighted age
    zu = input(' Compute flux weighted age in the galaxy rest frame at z [0] = ')
    if len(zu) <= 0:
        zu = 0.
    else:
        zu = float(zu)
    bt.zu = zu
    print()

def sfr(tp):
    # Returns SFR at time tp
    sfr=0.
    # Check for tcut
    if tp > tcut:
        return sfr
    # Check for wrong call
    if tp < 0:
        print(' SFR called with tp < 0, ' + str(tp))
        quit()
    # Select SFR according to io
    if io == 0:
        if tp == 0:
            sfr=1.
    elif io == 1:
        if tau > 0:
            sfr=(1. + epsilon*pgas(tp))*math.exp(-tp/tau)/tau
        elif tau < 0:
            sfr=(1./(1.-math.exp(1.)) + epsilon*pgas(tp))*math.exp(-tp/tau)/tau
    elif io == 2:
        if tp <= tau:
            sfr=1./tau
    elif io == 3:
        sfr=tau
    elif io == 4:
        if tau > 0:
            sfr = math.exp(-tp/tau)/tau
        elif tau < 0:
            # sfr = exp(-tp/tau)/tau/(1.-exp(1.))
            # Modified on June 6, 2017 to normalize to total mass = 1 at t = 13.6 Gyr
            # epsilon assumed to be = 0.
            tu = 13.6E9
            sfr = (math.exp(-tp/tau)-1.) / ( (1.-math.exp(-tu/tau))*tau - tu )
    elif io == 6:
        sfr=tp*math.exp(-tp/tau)/tau**2
    elif io == 7:
        sfr=usrsfr(tp,'')
    elif io == 8:
        sfr=2./tau*(1.-tp/tau)
        sfr=max(sfr,0.)
    elif io == 9:
        if ampr > 0:
            # If burst has been added to double exponential SFR, use
            sfr=usrsfr(tp,'')
        else:
            # Use analytic double exponential, joined at t = tauj
            if tp <= tauj:
                sfr = math.exp(-tp/tau1)/tau1
            else:
                sfr = con2*math.exp(-tp/tau2)
    return sfr

def sfrpar():
    print(' Choose SFR: 0 = SSP (Delta Burst = zero length burst)')
    print('             1 = Exponential (enter Tau)')
    print('            -1 = Exponential (enter mu_SFR parameter)')
    print('             2 = Single Burst of finite length')
    print('             3 = Constant')
    print('             4 = Delayed')
    print('             5 = Linearly decreasing')
    print('             6 = Read SFR(t) from ASCII file')
    print('             7 = Single or double exponential + burst (after Chen et al.)')
    io = input('    Choice = ')
    if io == '':
        io = '0'
    return int(io)

def csp_galaxev(file1):
    # python version of csp_galaxev code
    global f, t, iread, kf, fd, kp, s1, to, so
    global io, tau, ans, epsilon, tcut, tmax, ad, tv, mu, zu

    def CFredden(x,y,t,f,tv,mu):
        # Attenuate by dust if requested (use Charlot & Fall, 2000, ApJ, 539, 718)
        tauv = float(tv)
        d_mu = float(mu)
        p = -1
        if tauv <= 0:
            return y
        # Compute attenuation vs wavelength for young and old population components
        ay = []
        ao = []
        tauo = d_mu*tauv
        for i in range(len(x)):
            ay.append(fext(x[i],tauv))
            ao.append(fext(x[i],tauo))
        # Attenuate SED
        for j in range(len(t)):
            l = j+1
            for i in range(len(x)):
                if t[j] <= 1.E7:
                    y[i][l] = ay[i]*y[i][l]
                else:
                    y[i][l] = ao[i]*y[i][l]
            # Report percent done
            p = cn.percent(j,len(t),' Attenuating model SED, TauV = ' + str(tv) + ', mu = ' + str(mu) + ' --> ' + f,p)
        print()
        return y

    def fext(x,tv):
        # Attenuation function used in Charlot and Fall (2000) model
        if tv == 0:
            return 1
        tau=tv*( (5500./x)**0.7 )
        fext=math.exp(-tau)
        return fext

    def zflx():
        # Get fluxes needed to compute flux-weighted age
        global fd,ow
        kw = [120, 121, 122, 123, 124, 57]
        # Read filter filter file and build filter arrays at z = zu
        fd = fl.qfilters(x,kw,0.,zu,1)
        bt.ffd = fd
        # Compute flux in each band at all ages
        ip = len(fd[0])
        fz = []
        for k in range(len(t)):
            # Get sed corresponding to age a
            y = bc.ft(t[k],y1,t)
            v = []
            for i in range(ip):
                filtern = fl.qfilter_n(i,x,y,zu)
                v.append(filtern)
            fz.append(v)
        return fz

    def fzi(age,f,t):
        # Return flux in record corresponding to age t interpolating if necessary
        age = bc.age_yr(age)
        i1 = np.searchsorted(t,age)
        tlog = cn.tlog
        if t[i1] == age:
            # Record of desired age exists
            y1 = np.array(f[i1])
            ya = y1*t[i1]
            yl = y1*tlog(t[i1])
            yd = y1
            ma = t[i1]
            ml = tlog(t[i1])
            md = 1
        else:
            # Interpolate record
            i2=i1-1
            if age > 0 and t[i1] > 0 and t[i2] > 0:
                a1=np.log10(t[i2]/age)/np.log10(t[i2]/t[i1])
            else:
                a1=(t[i2]-age)/(t[i2]-t[i1])
            a2 = 1.-a1
            y1 = np.array(f[i1])
            y2 = np.array(f[i2])
            ya = a1*y1*t[i1]       + a2*y2*t[i2]
            yl = a1*y1*tlog(t[i1]) + a2*y2*tlog(t[i2])
            yd = a1*y1             + a2*y2
            ma = a1*t[i1]          + a2*t[i2]
            ml = a1*tlog(t[i1])    + a2*tlog(t[i2])
            md = a1 + a2
        ot = np.concatenate((ya, yl, yd, np.array([ma, ml, md])))
        return ot

    def convolve(ic,t):
        # Computes sed at age tb(k) by performing convolution integral of the SED for an SSP and the
        # chosen SFR using trapezoidal rule.
        # Follows CONVOLVE_NEW routine in fortran version (written by G. Bruzual on 25-March-1999)
        # The SFR is used as given (from t=0 to t=age) and the SSP spectra are interpolated at the required age.
        # This allows for arbitrarily short bursts to be described correctly.
        global to,so,wa

        def ts(ic):
            # Build time scale for accurate integration
            import numpy.ma as ma
            if ic == 0:         # This avoid zero divide in t = 0
                ic = 1          # Fortran version does not start from t = 0
            age  = t[ic]
            taux = []
            for k in range(ic):
                taux.append(t[k])
                taux.append(age-t[k])
            # Add more time steps under option 7 (user sfr entered in table)
        #   if io == 7:
        #       for i in range(len(time)):
        #           if age >= time[i]:
        #               taux.append(age-time[i])
        #       taux.append(time[i])
            # Add more steps if tcut < tmax
            if age > tcut and tcut < tmax:
                for k in range(ic):
                    if t[k] <= tcut:
                        last=k
                tlast=t[last]
                tnext=t[last+1]
                dn=0.0005*tcut
                for i in range(100000):
                    tlast=tlast+dn
                    if tlast <= tnext:
                        taux.append(age-tlast)
                        taux.append(tlast)
                    else:
                        break
            # Sort array taux and suppress duplicate entries, then reverse
            taux = np.unique(taux)
            # Suppress entries older than age
            taux = ma.masked_greater(taux, age)
            # Reverse array
            taux = np.flip(taux)
            return age,taux

        def dt(k):
            # Returns integration step dt
            if k == 0:
                dt=taux[0]-taux[1]
            elif k == kn-1:
                dt=taux[kn-2]-taux[kn-1]
            else:
                dt=taux[k-1]-taux[k+1]
            return dt

        # Construct integration time scale and perform convolution integral
        if io > 0:
            age, taux = ts(ic)
            kn = len(taux)
            y  = np.zeros(len(y1))
            b  = np.zeros(len(kp))
            wa = np.zeros(21)
            if kn==1:
                return y,b.tolist()
            for k in range(kn):
                sr=sfr(age-taux[k])
                if sr > 0:
                    tx=taux[k]
                    # Compute weight to assign to sed
                    ww=sr*dt(k)/2.
                    # Convolve SED
                    y = y + ww*bc.fi(tx,y1,t1)
                    # Convolve physical properties
                    bj = []
                    for j in range(len(kp)):
                        bj.append(ww*bc.pi(tx,p1,m1,t1,kp[j]))
                    b = b + np.array(bj)
                    # Convolve flux weighted age
                    wa = wa + ww*fzi(tx,fz,t1)

            # Store amount of processed gas
            if io == 1:
                xstr = b[10]
                xrm  = b[11]
                xstr = min(1,xstr)
                # unprocessed gas so far:
                ugas = math.exp(-age/tau)
                # processed gas = gas formed into stars - mass in stars - remnants
                prgas = 1. - ugas - xstr - xrm
                prgas = max(0,prgas)
                to.append(age)
                so.append(prgas)
                bt.to = to
                bt.so = so

        else:
            # SED
            y = ft(t1[ic],y1,t1)
            # Physical properties
            bj = []
            for j in range(len(kp)):
                bj.append(pt(t1[ic],p1,m1,t1,kp[j]))
            b = np.array(bj)
            # Flux weighted age
            wa = fzi(t1[ic],fz,t1)

        bt.wa = wa
        return y,b.tolist()

    # Init variables
    bt.iread = True	# Read filter file only on first call
    p=-1		# Needed by 'percent' function on first call
    to = [] ; so = []   # To store mass of recycled gas

    # Ask for input parameters to run csp_galaxev program
    if not jupy:
        bs.bchead()
        f1 = input(' BC_GALAXEV SSP sed in file [' +     bs.lnam(file1)  + '] = ')
        if len(f1) <= 0:
            f1 = bs.lnam(file1)
        else:
            f1 = bs.zrep(file1,f1)
    else:
            f1, fo = cspwidgtparam(file1)
            file1 = f1

    # Read BC/CB model in fits file
    x,y1,t1,e1,s1,m1,p1 = bc.bcnew(f1)
    bt.s1 = s1
    t=t1

    # Get SFR params
    if not jupy:
        cspparams()

        # Output file name
        fo = input(' Output file name = ')
        if fo=='q':
            quit()
        if fo == '' and bt.sn != '':
            fo = bs.lnam(f1).replace('.fits','').replace('ssp',bt.sn).strip()

    # Redden SED if requested
    if ad:
        y1 =CFredden(x,y1,t,f1,tv,mu)

    for i in range(len(t1)):
        if bt.iread:
            # Init output tables
            bt.kf,kp = cn.opent(0)
            # Read filter file and compute fluxes at z = zu
            fz = zflx()
            # Build filter arrays at z = 0
            bt.ffd = fl.qfilters(x,kf,0.,0,1)
            # Store model SED wavelength array
            bt.td1.add_column(np.array(x))
            # Do only on first call
            other = False

        # Perform convolution integral
        y,b = convolve(i,t1)

        # SFR at age t1[i]
        sf = sfr(t1[i])

        # Compute model rest-frame properties at this age
        cn.rf_props(t1[i],x,y,b,sf,0)

        # Store model SED = y
        bt.td1.add_column(np.array(y))

        # Report percent done
        p = cn.percent(i,len(t),' CSP_GALAXEV --> ' + fo,p)

    # Write results to fits file
    bc.wrfits(fo,0)

import os
import sys
import ipywidgets        as widgets
import numpy             as np
import builtins          as bt
import matplotlib.pyplot as plt
import bcfilt.bcfilt     as fl
import bcfits.bcfits     as bc
import basics.basics     as bs
import bchead.bchead     as hd
import cosmol.cosmol     as cl
import cspmod.cspmod     as csp
import addmod.addmod     as add
import bcevol.bcevol     as evl

def pyGALAXEVcl(a):
    # Read BC/CB models from fits table and call plotting functions
    global jupy, save
    jupy = False
    save = False
    bt.save = save
    bt.jupy = jupy
    bt.j  = False
    bt.ss = False
    bt.sc = False

    # Execute selected option - Command Line version
    n = pinit(a)
    while True:
        # Select plot option
        if n==0:
            # Ask for name of up to 2 fits files
            print()
           #n,bt.om = mfile()
            n,om[0],om[1] = mfile()
            bt.om = om
            n = rf(0,om)
        n, o, q = menu(0,n)
        fplot(n,o,q)

def sspFileName():
    # Builds file name corresponding to selected BC/CB model
  
    bc = {'None':['0.017','0','0.0001','0.0002','0.0005','0.001','0.002','0.004','0.006','0.008','0.010','0.014','0.017','0.02','0.03','0.04','0.06'],
          'BC19':['0.017','0','0.0001','0.0002','0.0005','0.001','0.002','0.004','0.006','0.008','0.010','0.014','0.017','0.02','0.03','0.04','0.06'],
          'CB19':['0.017','0','0.0001','0.0002','0.0005','0.001','0.002','0.004','0.006','0.008','0.010','0.014','0.017','0.02','0.03','0.04','0.06'],
         #'BC22':['0.017','0.0001','0.0005','0.001','0.002','0.004','0.006','0.008','0.010','0.014','0.017','0.02','0.03','0.04'],
         #'CB22':['0.017','0.0001','0.0005','0.001','0.002','0.004','0.006','0.008','0.010','0.014','0.017','0.02','0.03','0.04'],
          'BC03':['0.02','0.0001','0.0004','0.004','0.008','0.02','0.05'],
          'CB03':['0.019','0.0004','0.001','0.004','0.008','0.019','0.03'],
          'CB07':['0.02','0.0001','0.0004','0.004','0.008','0.02','0.05']}
    imf = ['Chabrier', 'Kroupa', 'Salpeter', 'Vaz_x0.3*', 'Vaz_x0.8*', 'Vaz_x1.0*', 'Vaz_x1.3*', 'Vaz_x1.5*', 'Vaz_x1.8*', 'Vaz_x2.0*', 'Vaz_x2.3*', 'Vaz_x2.8*', 'Vaz_x3.3*']
    mup = ['10*','100','300*','600*']
    atl = ['Miles','Miles+','IndoUS','Stelib','BaSeL']
    zW  = bc['CB19']
    dup = r'\(M_{UP}\)'

    # Define widget grid
    grid = widgets.GridspecLayout(2,9)

    # Minimal instructions
    ins = '<span style="color:blue">Use the dropdown menus to build the file name corresponding to the desired model(s). Values marked (*) are available only for the <b>CB19</b> models. The command <b>model ID</b> provides details on the model ingredients.</span>'
    grid[0,0:] = widgets.HTML(value=ins, placeholder='', description='',)

    # Dropdown menu's to build model file name
    grid[1,1] = widgets.Dropdown(options = bc.keys(), description='Model', value='CB19',layout={'width': 'max-content'})
    grid[1,2] = widgets.Dropdown(description='Z',layout={'width': 'max-content'})
    grid[1,3] = widgets.Dropdown(options=imf, description='IMF', value=imf[0], rows=len(imf), interactive=True, layout={'width': 'max-content'})      # Create dropdown widget to select IMF
    grid[1,4] = widgets.Dropdown(options=mup, description=dup, value=mup[1], rows=len(mup), interactive=True, layout={'width': 'max-content'})      # Create dropdown widget to select Mup
    grid[1,5] = widgets.Dropdown(options=atl, description='Library', value=atl[0], rows=len(atl), interactive=True, layout={'width': 'max-content'}) # Create dropdown widget to select Spectral Library
    mW  = grid[1,1]
    zW  = grid[1,2]
    imf = grid[1,3]
    mup = grid[1,4]
    lib = grid[1,5]
    def fn(Model, Z, IMF, Mup, Library):
        zW.options = bc[Model] # Here is the trick, i.e. update zW.options based on model, namely modelW.value.
        fw = namewidget(Model,zW.value,IMF,Mup,Library)
        # Textbox widget to select file 1
        grid[1,0] = widgets.HTML(value=fw, placeholder='', description='File:',)
    out = widgets.interactive_output(fn,{'Model':mW, 'Z':zW, 'IMF':imf, 'Mup':mup, 'Library':lib})
    display(grid,out)

def widget1(cv,d1,d2,ph,off,f_i):
    # Unique widget for all plots
    print()
    print()
    print()

    # Define widget grid
    grid = widgets.GridspecLayout(2,4)

    # Create textbox widget to select quantities to plot
    grid[0,0] = widgets.Text(value='', description=d2, placeholder=ph, disabled=False)

    # Create dropdown widget for preselected sets to plot
    grid[1,0] = widgets.Dropdown(options=cv, description=d1, value=cv[0], rows=len(cv), interactive=True)

    # Create textbox widget to select model 1
    grid[0,1] = widgets.Text(value='', placeholder='Model filename', description='Model 1', disabled=False)

    # Create textbox widget to select model 2
    grid[1,1] = widgets.Text(value='', placeholder='Model filename', description='Model 2', disabled=False)

    # Create radio button widget to select savefig option
    pv = ['Save fig(s) to png file(s)', 'Do not save fig(s)']
    grid[0,2] = widgets.RadioButtons(options=pv, value=pv[1], description='Save:',disabled=False)

    # Create radio button widget to select report file reading
    # rv = ['Report reading file', 'Do not report']
    # grid[0,3] = widgets.RadioButtons(options=rv, value=rv[1], description='Reading:',disabled=False)

    # Read widgets
    o0 = grid[0,0]		# Text menu
    o1 = grid[1,0]		# Dropdown menu
    o4 = grid[0,1]		# Model 1
    o5 = grid[1,1]		# Model 2
    o6 = grid[0,2]		# Save Fig option
   #o7 = grid[0,3]		# Report reading file (superseded)
    def f(o0,o1,o4,o5,o6):
        global ff,nn,rr,r1,r2,save
        ff = f_i
       #rr = rv.index(o7)
        ss = pv.index(o6)
        if ss==0:
            save=True
        else:
            save=False
        if len(o0) > 0:
            o0 = o0.replace(",", " ")
            bt.om[7] = o0.split()
            nn = off - 1
        else:
            nn = cv.index(o1) + off
        if len(o4) > 0:
            bt.om[0] = o4
            r1 = 1
        else:
            bt.om[0] = file1
            r1 = 0
        if len(o5) > 0:
            bt.om[1] = o5
            r2 = 1
        else:
            bt.om[1] = file2
            r2 = 0
        bt.sc = save
        bt.rr = 0
    out = widgets.interactive_output(f,{'o0':o0,'o1':o1,'o4':o4,'o5':o5,'o6':o6})

   # Create button widget to control plots
    button = widgets.Button(description='Plot')
    button.style.button_color = 'cyan'
    output = widgets.Output()
    grid[1,2] = button
    def f2(b):
        with output:
            if r1+r2 > 0:
                n = rf(rr,om)
            ff(nn)
    display(grid, output)
    button.on_click(f2)

def widget2(off,f_i):
    # Widget for sed plots
    print()
    print()

    # Define widget grid
    grid = widgets.GridspecLayout(2,6)

    # Create textbox widget to select quantities to plot
    d2 = 'Age of SEDs' ; ph = 'Enter age in Gyr'
    grid[0,0] = widgets.Text(value='', description=d2, placeholder=ph, disabled=False)

    # Create dropdown widget for preselected sets to plot
    cv = ['t < 100 Myr', 't < 1 Gyr', 't < 14 Gyr']   ; d1 = 'Preselected'
    grid[1,0] = widgets.Dropdown(options=cv, description=d1, value=cv[2], rows=len(cv), interactive=True)

    # Create textbox widget to select plot limits
    d2 = 'Plot limits' ; ph = 'xmin xmax log(ymin) log(ymax)'
    grid[0,1] = widgets.Text(value='', description=d2, placeholder=ph, disabled=False)

    # Create textbox widget to select normalization wavelength
    d2 = 'Normalize' ; ph = 'at wavelength (A)'
    grid[1,1] = widgets.Text(value='', description=d2, placeholder=ph, disabled=False)

    # Create textbox widget to select model 1
    ph = 'Model filename'
    ph = bs.lnam(file1)
    grid[0,2] = widgets.Text(value='', placeholder=ph, description='Model 1', disabled=False)

    # Create textbox widget to select model 2
    ph = 'Model filename'
    ph = bs.lnam(file2)
    grid[1,2] = widgets.Text(value='', placeholder=ph, description='Model 2', disabled=False)

    # Create radio button widget to select savefig option
    pv = ['Save fig(s) to png file(s)', 'Do not save fig(s)']
    grid[0,4] = widgets.RadioButtons(options=pv, value=pv[1], description='Save:',disabled=False)

    # Create radio button widget to select report reading file option (superseded)
    #rv = ['Report reading file', 'Do not report']
    #grid[0,4] = widgets.RadioButtons(options=rv, value=rv[1], description='Reading:',disabled=False)

    # Minimal instructions
    grd2 = widgets.GridspecLayout(5,6)
    v = '&nbsp;'
    def s(a,b):
        d = '<font size="2"><b>' + a + '</b><font size="2">( ' + b + ')' + 5*v + '<font size="2">'
        return d
    ins = '<b>Plot_SEDs</b>: to plot SEDs of age <b>t</b> in the wavelength range <b>[W1,W2]</b> follow the notation below to select the desired parameters and click the <b>Plot</b> button'
    grd2[0,0:]  = widgets.HTML(value=ins, placeholder='', description='')

    ins = 21*v + 'Age of SEDs:' + 3*v + s('-t','set of SEDs with age < t') + s('t1 t2 t3...','age of individual SEDs') + s('t1:t2','all SEDs in age range t1:t2') + \
          s('-N1:N2','all SEDs with record number in range N1:N2')
    grd2[1,0:]  = widgets.HTML(value=ins, placeholder='', description='')

    ins = 46*v + '<b>-t1;-t2</b>' + 8*v + '<b>t1;t2</b>' + 8*v + '<b>-t1;t2</b>' + 8*v + '<b>t1;-t2</b>' + 8*v + '<b>t1 t2 t3...;t4 t5 t6...</b>' + 8*v + \
          '<b>t1:t2;t3:t4</b>' + 8*v + '<b>-N1:N2;-N3:N4</b>' + 8*v + '(will plot different SEDs on each panel when plotting two models)'
    grd2[2,0:] = widgets.HTML(value=ins, placeholder='', description='')

    ins = 21*v + 'Plot limits:' + 8*v + '<b>xmin xmax</b> (leave empty for full range)' + 9*v + '<b>-xmin</b> to plot <b>Fnu</b>' + 9*v + '<b>log(ymin) log(ymax)</b> are optional'
    grd2[3,0:]  = widgets.HTML(value=ins, placeholder='', description='')
    display(grd2)

    # Read widgets
    o0 = grid[0,0]		# Age array
    o1 = grid[1,0]		# Dropdown menu
    o2 = grid[0,1]		# Plot limits
    o3 = grid[1,1]		# Normalization wavelength
    o4 = grid[0,2]		# Model 1
    o5 = grid[1,2]		# Model 2
    o6 = grid[0,4]		# Save Fig option
   #o7 = grid[0,4]		# Report reading file (superseded)

    def f1(o0,o1,o2,o3,o4,o5,o6):
        global nn,ff,rr,r1,r2,save
        ff = f_i
        ss = pv.index(o6)
       #rr = rv.index(o7)
        if ss==0:
            save=True
        else:
            save=False
        if len(o0) > 0:
            o0 = o0.replace(',', '|')
            o0 = o0.replace(' ', '|')
            om[7] = o0
            nn = off
        else:
            nn = cv.index(o1) + off + 1
        o2    = o2.replace(',', ' ')
        om[8] = o2.split()
        om[9] = o3
        if len(o4) > 0:
            om[0] = o4
            r1 = 1
        else:
            o4 = file1
            om[0] = o4
            r1 = 1
        if len(o5) > 0:
            om[1] = o5
            r2 = 1
            if nn < 4:
                nn = nn+4
        else:
            o5 = file2
            om[1] = o5
            r2 = 1
        bt.om = om
        bt.r1 = r1
        bt.r2 = r2
        bt.nn = nn
        bt.ss = save
        bt.rr = 0
    out = widgets.interactive_output(f1,{'o0':o0,'o1':o1,'o2':o2,'o3':o3,'o4':o4,'o5':o5,'o6':o6})

   # Create button widget to print time steps
   #button= widgets.Button(description='List time steps',layout=widgets.Layout(width='90%', height='30px'))
    button= widgets.Button(description='List time steps')
    button.style.button_color = 'yellow'
    output = widgets.Output()
    grid[1,4] = button
    def g2(b):
        with output:
            oo = [om[0],'s']
            bc.bcf2t(oo)
    button.on_click(g2)

   # Create button widget to control plots
   #button = widgets.Button(description='Plot',layout=widgets.Layout(width='90%', height='30px'))
    button = widgets.Button(description='Plot')
    button.style.button_color = 'cyan'
    output = widgets.Output()
    grid[1,3] = button
    def f2(b):
        with output:
            if r1+r2 > 0:
                n = rf(rr,om)
            ff(nn)
    display(grid, output)
    button.on_click(f2)

def widgetrfphot():
    # Create widget to run rf_phot program from Jupyter notebook

    # Define widget grid
    grid = widgets.GridspecLayout(6,9)

    # Copyright + revision date
   #grid[0,0:] = widgets.HTML(value=cry, placeholder='', description='')

    # Minimal instructions
    ins = '<b>rf_phot</b>: compute photometric magnitudes in the galaxy rest frame in selected filter bands'
    grid[1,0:] = widgets.HTML(value=ins, placeholder='', description='')
    v = '&nbsp;'
    ins = 16*v + 'select desired <b>galaxy age, redshift, extinction Av</b>, and <b>filter set</b> (or leave empty for default values) and click the <b>run rf_phot</b> button'
    grid[2,0:] = widgets.HTML(value=ins, placeholder='', description='')

    # BC model
    file1 = bt.fw
    dn = 'Input model' ; ph = bs.lnam(file1) 	# ph = 'Enter input file name'
    grid[3,0] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Output file (enter output file name)
    dn = 'Output File' ; ph = 'Enter output file name'
    grid[4,0] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Create textbox widget to select age of SEDs
    d2 = 'Age of SSP' ; ph = 'Enter age in Gyr [all]'
    grid[3,1] = widgets.Text(value='', description=d2, placeholder=ph, disabled=False)

    # galaxy redshift
    dn = 'z' ; ph = 'Enter galaxy redshift [0]'
    grid[4,1] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # extinction Av
    dn = 'Av' ; ph = 'Enter Av [0]'
    grid[5,1] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Select filters by number
    dn = 'Filters' ; ph = 'Filters'
    grid[3,2] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Select magnitude system
    cv = ['VEGAmag', 'ABmag', 'STmag'] ; d1 = 'Magnitude'
    grid[4,2] = widgets.Dropdown(options=cv, description=d1, value=cv[1], rows=len(cv), interactive=True)

    # Create dropdown widget for preselected filter sets
    cv = ['Johnson (Vega mag)', 'SDSS (AB mag)', 'HST/JWST (ST mag)', 'J-PLUS (AB mag)', 'JPAS (AB mag)' ]   ; d1 = 'Preselected'
    grid[5,2] = widgets.Dropdown(options=cv, description=d1, value=cv[1], rows=len(cv), interactive=True)

    # Read widgets
    o0 = grid[3,0]		# Input file
    o1 = grid[4,0]		# Output file
    o2 = grid[3,2]		# Filter set (entered by user)
    o3 = grid[4,2]		# Magnitude system (dropdown menu)
    o4 = grid[5,2]		# Preselected filter set (dropdown menu)
    o5 = grid[3,1]		# Galaxy age
    o6 = grid[4,1]		# Galaxy redshift
    o7 = grid[5,1]		# Extinction Av

    def f1(o0,o1,o2,o3,o4,o5,o6,o7):
        global r1,r2
        rm[0] = o0
        rm[1] = o1
        rm[2] = o2
        rm[3] = o3
        rm[4] = o4
        rm[5] = o5
        rm[6] = o6
        rm[7] = o7
        # Input model
        if len(o0) > 0:
            file = o0
        else:
            file = file1
        # Output file
        if len(o1) > 0:
            of = '-o ' + o1
        else:
            of = '-o'
        # Age, z and Av
        if len(o5) > 0:
            tx = 't' + o5
        else:
            tx = ''
        if len(o6) > 0:
            zx = 'z' + o6
        else:
            zx = ''
        if len(o7) > 0:
            av = 'av' + o7
        else:
            av = ''
        # Preselected filter set
        if 'Johnson' in o4:
            fds = 'john'
        elif 'SDSS' in o4:
            fds = 'sdss'
        elif 'HST' in o4:
            fds = 'stmag'
        elif 'PLUS' in o4:
            fds = 'jplus'
        elif 'JPAS' in o4:
            fds = 'jpas'
        else:
            fds = 'sdss'
        # User selected filters
        if len(o2) > 0:
            o2  = o2.replace(',',' ')
            o2  = 'user ' + o3 + ' ' + o2
            fds= o2.split()
        bt.rf = [ file, fds, tx, zx, av, of ]
    out = widgets.interactive_output(f1,{'o0':o0,'o1':o1,'o2':o2,'o3':o3,'o4':o4,'o5':o5,'o6':o6,'o7':o7})

   # Create button widget to control plots
    button = widgets.Button(description='run rf_phot')
    button.style.button_color = 'lightblue'
    output = widgets.Output()
    grid[5,4] = button
    def f2(b):
        with output:
            if len(rm[0]) > 0:
                bt.om[0] = rm[0]
                n = cf(bt.om)
                r1 = 1
            else:
                bt.om[0] = file1
                r1 = 0
            fplot(2,'5',1)
    display(grid, output)
    button.on_click(f2)

def widgetofphot():
    # Create widget to run of_phot program from Jupyter notebook

    # Define widget grid
    grid = widgets.GridspecLayout(8,9)

    # Copyright + revision date
   #grid[0,0:] = widgets.HTML(value=cry, placeholder='', description='')

    # Minimal instructions
    ins = '<b>of_phot</b>: compute photometric magnitudes in the observer frame in selected filter bands'
    grid[1,0:] = widgets.HTML(value=ins, placeholder='', description='')
    v = '&nbsp;'
    ins = 16*v + 'select desired <b>galaxy age, filter number(s)</b>, and <b>cosmological model parameters</b> (leave empty for default values) and click the <b>run of_phot</b> button'
    grid[2,0:] = widgets.HTML(value=ins, placeholder='', description='')

    # BC model
    file1 = bt.fw
    dn = 'Input model' ; ph = bs.lnam(file1)			# 'GALAXEV model file name'
    grid[3,0] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Output file (enter output file name)
    dn = 'Output File' ; ph = 'Enter output file name'
    grid[4,0] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Hubble constant
    dn = r'\(H_0\)' ; ph = '71 km/s/Mpc'
    grid[5,0] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Omega matter
    dn = r'\(\Omega_m\)' ; ph = '0.27'
    grid[6,0] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Omega_lambda
    dn = r'\(\Omega_\Lambda\)' ; ph = '0.73'
    grid[7,0] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # extinction Av
    dn = 'Av' ; ph = 'Enter Av [0]'
    grid[4,1] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Age of galaxies
    dn = r'\(t_g\)' ; ph = '13.5 Gyr (age of galaxy)'
    grid[5,1] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Select filters by number
    dn = 'Filters' ; ph = 'Filters'
    grid[5,2] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Select magnitude system
    cv = ['VEGAmag', 'ABmag', 'STmag'] ; d1 = 'Magnitude'
    grid[6,2] = widgets.Dropdown(options=cv, description=d1, value=cv[1], rows=len(cv), interactive=True)

    # Create dropdown widget for preselected filter sets
    cv = ['Johnson (Vega mag)', 'SDSS (AB mag)', 'HST/JWST (ST mag)', 'J-PLUS (AB mag)', 'JPAS (AB mag)' ]   ; d1 = 'Preselected'
    grid[7,2] = widgets.Dropdown(options=cv, description=d1, value=cv[1], rows=len(cv), interactive=True)

    # Read widgets
    o0 = grid[3,0]		# Input file name
    o1 = grid[4,0]		# Output file name
    o2 = grid[5,0]		# Ho
    o3 = grid[6,0]		# Omega_m
    o4 = grid[7,0]		# Omega_lambda
    o5 = grid[5,2]		# Filter set (entered by user)
    o6 = grid[6,2]		# Magnitude system (dropdown menu)
    o7 = grid[7,2]		# Preselected filter set (dropdown menu)
    o8 = grid[4,1]		# Extinction Av
    o9 = grid[5,1]		# Galaxy age

    def f(o0,o1,o2,o3,o4,o5,o6,o7,o8,o9):
        global r1,tu,zf
        rm[1] = o1
        rm[2] = o2
        rm[3] = o3
        rm[4] = o4
        rm[5] = o5
        rm[6] = o6
        rm[7] = o7
        rm[8] = o8
        rm[9] = o9
        if len(o0) > 0:
            file = o0
            r1 = 1
        else:
            file = file1
            r1 = 0
        # Define default value of cosmological parameters
        h, q, clambda, tu, ttg, zf, omega, omega_lambda = bt.cosmol
        if len(rm[2]) > 0:
            h = float(rm[2])
        if len(rm[3]) > 0:
            omega = float(rm[3])
        if len(rm[4]) > 0:
            omega_lambda = float(rm[4])
        if len(rm[9]) > 0:
            ttg = float(rm[9])
        clambda,q = cl.cosmol_c(h,omega,omega_lambda)

        # z of galaxy formation
        zf = cl.zx(ttg,h,q,clambda)
        vl = str('{:.1f}'.format(zf))
        ph = vl + ' (z of galaxy formation)'
        dn = r'\(z_f\)'
        grid[6,1] = widgets.Text(value='', description=dn, placeholder=ph, disabled=True)
        bt.cosmol2 = h, q, clambda, tu, ttg, zf, omega, omega_lambda

        # age of universe
        tu = cl.tuniverse(h,q,0.,clambda)
        vl = str('{:.1f}'.format(tu))
        ph = vl + ' Gyr (age of universe)'
        dn = r'\(t_u\)'
        grid[7,1] = widgets.Text(value='', description=dn, placeholder=ph, disabled=True)

        # Filter set
        if 'Johnson' in o7:
            fds = 'john'
        elif 'SDSS' in o7:
            fds = 'sdss'
        elif 'HST' in o7:
            fds = 'stmag'
        elif 'PLUS' in o7:
            fds = 'jplus'
        elif 'JPAS' in o7:
            fds = 'jpas'
        else:
            fds = 'sdss'
        # User selected filters
        if len(o5) > 0:
            o5  = o5.replace(',',' ')
            o5  = 'user ' + o6 + ' ' + o5
            fds= o5.split()
        if len(o8) > 0:
            av = 'av' + o8
        else:
            av = ''
        if len(o1) > 0:
            bt.of = [ file, fds, av, '-o', o1 ]
        else:
            bt.of = [ file, fds, av ]

    out = widgets.interactive_output(f,{'o0':o0,'o1':o1,'o2':o2,'o3':o3,'o4':o4,'o5':o5,'o6':o6,'o7':o7,'o8':o8,'o9':o9})

   # Create button widget to control flow
    button = widgets.Button(description='run of_phot')
    button.style.button_color = 'lightblue'
    output = widgets.Output()
    grid[7,4] = button
    def f2(b):
        with output:
            if len(rm[0]) > 0:
                bt.om[0] = rm[0]
                n = cf(bt.om)
                r1 = 1
            else:
                bt.om[0] = file1
                r1 = 0
            fplot(2,'6',1)
    display(grid, output)
    button.on_click(f2)

def widgetbcf2t():
    # Create widget to run bcfits2txt program from Jupyter notebook

    # Define widget grid
    grid = widgets.GridspecLayout(9,9)

    # Copyright + revision date
   #grid[0,0:] = widgets.HTML(value=cry, placeholder='', description='')

    # Minimal instructions
    ins = '<b>bcfits2txt</b>: create ascii files with selected model records or new files in various legacy formats'
    grid[1,0:] = widgets.HTML(value=ins, placeholder='', description='')
    b = '&nbsp;'
    ins = 19*b + ' select <b>records</b> to extract or <b>files</b> to create following the indicated notation and click the <b>run bcfits2txt</b> button'
    grid[2,0:] = widgets.HTML(value=ins, placeholder='', description='')

    # BC model
    dn = 'Input model' ; ph = bs.lnam(file1) 	# ph = 'Enter input file name'
    grid[3,0] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Records to extract
    dn = 'Records' ; ph = 't1 t2 t3...'
    grid[4,0] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Files to create
    dn = 'Files' ; ph = 't'
    grid[4,2] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # instructions
    v = '&nbsp;'
    ins1 = 7*v + '<b>t1 t2 t3...</b>' + 5*v + '(age of SEDs in Gyr)' + 75*v + '<b>t</b>' + 2*v + '(write tables with different properties of SSP model)'
    ins2 = 6*v + '<b>-N1 N2 N3...</b> (record number of SEDs. Note - sign in front of N1)' + 23*v + '<b>i</b>' + 2*v + '(create .ised fortran readable binary file compatible with BC03 model release)'
    ins3 = 7*v + '<b>N1:N2</b>' + 11*v + '(all SEDs with N1 <= N <= N2)' + 58*v + '<b>c</b>' + 2*v + '(create .ised_ASCII text file compatible with BC03 model release)'
    ins4 = 7*v + '<b>a</b>' + 20*v + '(all SEDs in fits file)'
    grid[5,0:] = widgets.HTML(value=ins1, placeholder='', description='')
    grid[6,0:] = widgets.HTML(value=ins2, placeholder='', description='')
    grid[7,0:] = widgets.HTML(value=ins3, placeholder='', description='')
    grid[8,0:] = widgets.HTML(value=ins4, placeholder='', description='')

    # Read widgets
    o0 = grid[3,0]		# Input file
    o1 = grid[4,0]		# Records to extract
    o2 = grid[4,2]		# Table selection
    if o2 == 'I':
        o2 = 'i'

    def f1(o0,o1,o2):
        global r1,r2
        rm[0] = o0
        rm[1] = o1
        rm[2] = o2
        if len(o0) > 0:
            bt.om[0] = o0
            r1 = 1
        else:
            bt.om[0] = file1
            r1 = 0
        if len(o1) > 0:
            bt.om[1] = o1
        if len(o2) > 0:
            bt.om[1] = o2
        bt.om = [bt.om[0], bt.om[1]]
    out = widgets.interactive_output(f1,{'o0':o0,'o1':o1,'o2':o2})

   # Create button widget to control plots
    button = widgets.Button(description='run bcfits2txt')
    button.style.button_color = 'lightblue'
    output = widgets.Output()
    grid[4,4] = button
    def f2(b):
        with output:
           #if len(rm[0]) > 0:
           #    bt.om[0] = rm[0]
           #    n = cf(bt.om)
           #    r1 = 1
           #else:
           #    bt.om[0] = file1
           #    r1 = 0
            fplot(2,'7',1)
    display(grid, output)
    button.on_click(f2)

def widgetadd():
    # Create widget to run add_bursts program from Jupyter notebook

    # Define widget grid
    grid = widgets.GridspecLayout(6,10)

    # Copyright + revision date
   #grid[0,0:] = widgets.HTML(value=cry, placeholder='', description='')

    # Minimal instructions
    ins = '<b>add_bursts</b>: compute the combined sed for 2 bursts of star formation. The sed corresponding to each burst must be computed first. Each burst is characterized by the beginning time (in Gyr)'
    grid[1,0:] = widgets.HTML(value=ins, placeholder='', description='')
    b = '&nbsp;'
    ins = 23*b + 'and the burst strength. Each individual burst can follow an arbitrary star formation law. Select burst parameters and click the <b>run add_burst</b> button'
    grid[2,0:] = widgets.HTML(value=ins, placeholder='', description='')

    # First model
    dn = '1st burst' ; ph = bs.lnam(file1) 		# 'Enter GALAXEV model file name'
    grid[3,0] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Second model
    dn = '2nd burst' 					# ; ph = 'Enter GALAXEV model file name'
    grid[4,0] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Output file (enter output file name)
    dn = 'Output File' ; ph = 'Enter output file name'
    grid[5,0] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Burst 1
    dn = 'Burst 1' ; ph = 'Beginning time (Gyr), strength'
    grid[3,1] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Burst 1
    dn = 'Burst 2' ; ph = 'Beginning time (Gyr), strength'
    grid[4,1] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Read widgets
    o0  = grid[3,0]		# First model
    o1  = grid[4,0]		# Second model
    o2  = grid[3,1]		# Burst 1 parameters
    o3  = grid[4,1]		# Burst 2 parameters
    o4  = grid[5,0]		# Output file name

    def f(o0,o1,o2,o3,o4):
        global r1,r2
        if len(o0) > 0:
            bt.om[0] = o0
            r1 = 1
        else:
            bt.om[0] = file1
            r1 = 0
        if len(o1) > 0:
            bt.om[1] = o1
            r2 = 1
        else:
            bt.om[1] = file1
            r2 = 0
        if len(o2) > 0:
            o2 = o2.replace(" ", ",")
            rm[0] = o2
        if len(o3) > 0:
            o3 = o3.replace(" ", ",")
            rm[1] = o3
        rm[2] = o4
    out = widgets.interactive_output(f,{'o0':o0,'o1':o1,'o2':o2,'o3':o3,'o4':o4})

   # Create button widget to control flow
    button = widgets.Button(description='run add_burst')
    button.style.button_color = 'lightblue'
    output = widgets.Output()
    grid[5,3] = button
    def f2(b):
        with output:
            if r1+r2> 0:
                n = rf(0,om)
            fplot(2,'1',1)
    display(grid, output)
    button.on_click(f2)

def widgetcmev():
    # Create widget to run cm_evolution program from Jupyter notebook

    # Define widget grid
    grid = widgets.GridspecLayout(7,9)

    # Copyright + revision date
   #grid[0,0:] = widgets.HTML(value=cry, placeholder='', description='')

    # Minimal instructions
    ins = '<b>cm_evolution</b>: compute photometric magnitudes and k-corrections in the observer frame in up to two selected filter bands as a function of redshift'
    grid[1,0:] = widgets.HTML(value=ins, placeholder='', description='')
    b = '&nbsp;'
    ins = 26*b + 'select desired <b>galaxy age, filter number(s)</b>, and <b>cosmological model parameters</b> (leave empty for default values) and click the <b>run cm_evolution</b> button'
    grid[2,0:] = widgets.HTML(value=ins, placeholder='', description='')

    # BC model
    dn = 'Model' ; ph = bs.lnam(file1)			# 'GALAXEV model file name'
    grid[3,0] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Hubble constant
    dn = r'\(H_0\)' ; ph = '71 km/s/Mpc'
    grid[4,0] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Omega matter
    dn = r'\(\Omega_m\)' ; ph = '0.27'
    grid[5,0] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Omega_lambda
    dn = r'\(\Omega_\Lambda\)' ; ph = '0.73'
    grid[6,0] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Age of galaxies
    dn = 'tg' ; ph = '13.5 Gyr (age of galaxy)'
    dn = r'\(t_g\)'
    grid[4,1] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # 1st filter
    dn = 'Filter 1' ; ph = '15'
    grid[5,2] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # 2nd filter
    dn = 'Filter 2' ; ph = 'Enter filter number'
    grid[6,2] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Read widgets
    o0  = grid[3,0]		# File name
    o1  = grid[5,2]		# Filter No 1
    o2  = grid[6,2]		# Filter No 2
    o3  = grid[4,1]		# Galaxy age
    o4  = grid[4,0]		# Ho
    o5  = grid[5,0]		# Omega_m
    o6  = grid[6,0]		# Omega_lambda

    def f(o0,o1,o2,o3,o4,o5,o6):
        global r1,tu,zf
        if len(o0) > 0:
            bt.om[0] = o0
            r1 = 1
        else:
            bt.om[0] = file1
            r1 = 0
        rm[1] = o1
        rm[2] = o2
        rm[3] = o3
        rm[4] = o4
        rm[5] = o5
        rm[6] = o6
        # Define default value of cosmological parameters
        h, q, clambda, tu, ttg, zf, omega, omega_lambda = bt.cosmol
        if len(rm[4]) > 0:
            h = float(rm[4])
        if len(rm[5]) > 0:
            omega = float(rm[5])
        if len(rm[6]) > 0:
            omega_lambda = float(rm[6])
        if len(rm[3]) > 0:
            ttg = float(rm[3])
        clambda,q = cl.cosmol_c(h,omega,omega_lambda)

        # age of universe
        tu = cl.tuniverse(h,q,0.,clambda)
        vl = str('{:.1f}'.format(tu))
        ph = vl + ' Gyr (age of universe)'
        dn = r'\(t_u\)'
        grid[5,1] = widgets.Text(value='', description=dn, placeholder=ph, disabled=True)

        # z of galaxy formation
        zf = cl.zx(ttg,h,q,clambda)
        vl = str('{:.1f}'.format(zf))
        ph = vl + ' (z of galaxy formation)'
        dn = r'\(z_f\)'
        grid[6,1] = widgets.Text(value='', description=dn, placeholder=ph, disabled=True)
        bt.cosmol2 = h, q, clambda, tu, ttg, zf, omega, omega_lambda

    out = widgets.interactive_output(f,{'o0':o0,'o1':o1,'o2':o2,'o3':o3,'o4':o4,'o5':o5,'o6':o6})

   # Create button widget to control flow
    button = widgets.Button(description='run cm_evolution')
    button.style.button_color = 'lightblue'
    output = widgets.Output()
    grid[6,4] = button
    def f2(b):
        with output:
            if r1> 0:
                n = cf(bt.om)
            fplot(2,'2',1)
    display(grid, output)
    button.on_click(f2)

def widgetcsp():
    # Create widget to run csp_galaxev program from Jupyter notebook

    # Define widget grid
    grid = widgets.GridspecLayout(9,3)

    # Copyright + revision date
   #grid[0,0:] = widgets.HTML(value=cry, placeholder='', description='')

    # Minimal instructions
    ins = '<b>csp_galaxev</b>: compute the properties of a composite stellar propulation that forms stars according to a given star formation law'
    grid[1,0:] = widgets.HTML(value=ins, placeholder='', description='')
    v = '&nbsp;'
    ins = 25*v + 'select desired <b>SFR</b> and othe model properties (leave empty for default values) and click the <b>run csp_galaxev</b> button'
    grid[2,0:] = widgets.HTML(value=ins, placeholder='', description='')

    # Output file (enter output file name)
    dn = 'Input SSP' ; ph = bs.lnam(file1) 	# ph = 'Enter input file name'
    grid[3,0] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Exponential SFR (enter tau)
    dn = 'Tau-model' ; ph = 'Enter e-folding time TAU (Gyr)'
    grid[3,1] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # c-model (enter burst duration)
    dn = 'c-model' ; ph = 'Enter duration of burst (Gyr)'
    grid[4,1] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Exponential SFR (enter tau)
    dn = 'Constant' ; ph = 'Enter SFR in Mo/yr'
    grid[5,1] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Delayed SFR (enter tau)
    dn = 'Delayed' ; ph = 'Maximum in SFR at time TAU (Gyr)'
    grid[6,1] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Linear SFR (enter tau)
    dn = 'Linear' ; ph = 'SFR = 0 at time TAU (Gyr)'
    grid[7,1] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Numerical SFR (enter table name)
    dn = 'Numerical' ; ph = 'Enter name of file with SFR'
    grid[8,1] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Dust model (enter tau_V and mu)
    dn = 'Dust model' ; ph = 'Enter tau_V, mu (Charlot & Fall)'
    grid[3,2] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Compute flux weighted age in the galaxy rest frame at z
    dn = 'WeightedAge' ; ph = 'Enter z of galaxy rest frame'
    grid[4,2] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Compute flux weighted age in the galaxy rest frame at z
    dn = 'Quench' ; ph = 'Make SFR = 0 at TCUT (Gyr)'
    grid[5,2] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Output file (enter output file name)
    dn = 'Output File' ; ph = 'Enter output file name'
    grid[4,0] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Read widgets
    o0  = grid[3,0]		# input file
    o1  = grid[3,1]		# tau model
    o2  = grid[4,1]		# c-model
    o3  = grid[5,1]		# constant sfr model
    o4  = grid[6,1]		# delayed sfr model
    o5  = grid[7,1]		# linear sfr model
    o6  = grid[8,1]		# numerical sfr model
    o7  = grid[3,2]		# dust model
    o8  = grid[4,2]		# flux weighted age
    o9  = grid[5,2]		# quench age
    o10 = grid[4,0]		# output file name

    def f(o0,o1,o2,o3,o4,o5,o6,o7,o8,o9,o10):
        global r1
        if len(o0) > 0:
            bt.om[0] = o0
            r1 = 1
        else:
            bt.om[0] = file1
            r1 = 0
        if len(o1) > 0:
            rm[1] = '1'
            rm[2] = o1
        elif len(o2) > 0:
            rm[1] = '2'
            rm[2] = o2
        elif len(o3) > 0:
            rm[1] = '3'
            rm[2] = o3
        elif len(o4) > 0:
            rm[1] = '4'
            rm[2] = o4
        elif len(o5) > 0:
            rm[1] = '5'
            rm[2] = o5
        elif len(o6) > 0:
            rm[1] = '6'
            rm[2] = o6
        if len(o7) > 0:
            o7 = o7.replace(",", " ")
            rm[3] = o7.split()
        rm[4] = o8
        rm[5] = o9
        rm[6] = o10
    out = widgets.interactive_output(f,{'o0':o0,'o1':o1,'o2':o2,'o3':o3,'o4':o4,'o5':o5,'o6':o6,'o7':o7,'o8':o8,'o9':o9,'o10':o10})

   # Create button widget to control flow
    button = widgets.Button(description='run csp_galaxev')
    button.style.button_color = 'lightblue'
    output = widgets.Output()
    grid[8,2] = button
    def f2(b):
        with output:
            if r1 > 0:
                n = rf(0,om)
            fplot(2,'0',1)
    display(grid, output)
    button.on_click(f2)

def widgetgpl():
    # Create widget to run galaxev_pl program from Jupyter notebook

    # Define widget grid
    grid = widgets.GridspecLayout(7,3)

    # Copyright + revision date
   #grid[0,0:] = widgets.HTML(value=cry, placeholder='', description='')

    # Minimal instructions
    ins = '<b>galaxevpl</b>: create ascii files with selected model records'
    grid[1,0:] = widgets.HTML(value=ins, placeholder='', description='')
    v = '&nbsp;'
    ins = 20*v + 'enter <b>age</b> of SEDs and desired <b>wavelength range</b> [W1,W2] in A (leave empty for default values) and click the <b>run galaxevpl</b> button'
    grid[2,0:] = widgets.HTML(value=ins, placeholder='', description='')

    # Create textbox widget to select model 1
    dn = 'Input model' ; ph = bs.lnam(file1) 	# ph = 'Enter input file name'
    grid[3,0] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Output file
    d2 = 'Output file' ; ph = 'Enter output file name'
    grid[3,1] = widgets.Text(value='', description=d2, placeholder=ph, disabled=False)

    # Create textbox widget to select age of SEDs
    d2 = 'Age of SEDs' ; ph = 'Enter ages in Gyr or -all'
    grid[4,0] = widgets.Text(value='', description=d2, placeholder=ph, disabled=False)

    # Create textbox widget to select wavelength range
    d2 = 'Range (A)' ; ph = 'W1, W2, Wo, Fo'
    grid[4,1] = widgets.Text(value='', description=d2, placeholder=ph, disabled=False)

    # Instructions
    ins = 148*v + 'Enter Wo,Fo for SEDs scaled to flux = Fo at lambda = Wo (default: no scaling)'
    grid[5,0:] = widgets.HTML(value=ins, placeholder='', description='')
    ins = 148*v + 'Enter W1 with a minus sign for SEDs listed as Fnu vs. lambda'
    grid[6,0:] = widgets.HTML(value=ins, placeholder='', description='')

    # Read widgets
    o0 = grid[3,0]		# Model file name
    o1 = grid[4,0]		# Age of SED'
    o2 = grid[4,1]		# Wavelength range
    o3 = grid[3,1]		# Output file name

    def f(o0,o1,o2,o3):
        rm[0] = o0
        rm[1] = o1
        rm[2] = o2
        rm[3] = o3
        bt.rm = rm
    out = widgets.interactive_output(f,{'o0':o0,'o1':o1,'o2':o2,'o3':o3})

   # Create button widget to control flow
    button = widgets.Button(description='run galaxevpl')
    button.style.button_color = 'lightblue'
    output = widgets.Output()
    grid[4,2] = button
    def f2(b):
        global om
        with output:
            if len(rm[0]) > 0:
                om[0] = rm[0]
                n = cf(om)
                r1 = 1
            else:
                om[0] = file1
                r1 = 0
            om[1] = ''
            bt.om = om
            bt.nu = False
            fplot(2,'4',1)
    display(grid, output)
    button.on_click(f2)

def widgetzmag():
    # Create widget to run zmag program from Jupyter notebook

    # Define widget grid
    grid = widgets.GridspecLayout(7,9)

    # Copyright + revision date
   #grid[0,0:] = widgets.HTML(value=cry, placeholder='', description='')

    # Minimal instructions
    ins = '<b>zmag</b>: compute photometric magnitudes and k-corrections in the observer frame in the selected filter band and redshift'
    grid[1,0:] = widgets.HTML(value=ins, placeholder='', description='')
    b = '&nbsp;'
    ins = 13*b + 'select desired <b>galaxy age, redshift, filter number</b>, and <b>cosmological model parameters</b> (leave empty for default values) and click the <b>run zmag</b> button'
    grid[2,0:] = widgets.HTML(value=ins, placeholder='', description='')

    # BC model
    dn = 'Model' ; ph = bs.lnam(file1)			# 'GALAXEV model file name'
    grid[3,0] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Hubble constant
    dn = r'\(H_0\)' ; ph = '71 km/s/Mpc'
    grid[4,0] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Omega matter
    dn = r'\(\Omega_m\)' ; ph = '0.27'
    grid[5,0] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Omega_lambda
    dn = r'\(\Omega_\Lambda\)' ; ph = '0.73'
    grid[6,0] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Age of galaxies
    dn = r'\(t_g\)' ; ph = '13.5 Gyr (age of galaxy today)'
    grid[4,1] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Filter
    dn = 'Filter' ; ph = '121'
    grid[5,2] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # galaxy redshift
    dn = 'z' ; ph = 'Enter galaxy redshift'
    grid[6,2] = widgets.Text(value='', description=dn, placeholder=ph, disabled=False)

    # Read widgets
    o0  = grid[3,0]		# File name
    o1  = grid[5,2]		# Filter
    o2  = grid[6,2]		# Galaxy redshift
    o3  = grid[4,1]		# Galaxy age
    o4  = grid[4,0]		# Ho
    o5  = grid[5,0]		# Omega_m
    o6  = grid[6,0]		# Omega_lambda

    def f(o0,o1,o2,o3,o4,o5,o6):
        global r1,tu,zf
        if len(o0) > 0:
            bt.om[0] = o0
            r1 = 1
        else:
            bt.om[0] = file1
            r1 = 0
        rm[1] = o1
        rm[2] = o2
        rm[3] = o3
        rm[4] = o4
        rm[5] = o5
        rm[6] = o6
        # Define default value of cosmological parameters
        h, q, clambda, tu, ttg, zf, omega, omega_lambda = bt.cosmol
        if len(rm[4]) > 0:
            h = float(rm[4])
        else:
            rm[4] = str(h)
        if len(rm[5]) > 0:
            omega = float(rm[5])
        else:
            rm[5] = str(omega)
        if len(rm[6]) > 0:
            omega_lambda = float(rm[6])
        else:
            rm[6] = str(omega_lambda)
        if len(rm[3]) > 0:
            ttg = float(rm[3])
        else:
            rm[3] = str(ttg)
        if len(rm[1]) == 0:
            rm[1] = '121'
        clambda,q = cl.cosmol_c(h,omega,omega_lambda)

        # age of universe
        tu = cl.tuniverse(h,q,0.,clambda)
        vl = str('{:.1f}'.format(tu))
        ph = vl + ' Gyr (age of universe)'
        dn = r'\(t_u\)'
        grid[5,1] = widgets.Text(value='', description=dn, placeholder=ph, disabled=True)

        # z of galaxy formation
        zf = cl.zx(ttg,h,q,clambda)
        vl = str('{:.1f}'.format(zf))
        ph = vl + ' (z of galaxy formation)'
        dn = r'\(z_f\)'
        grid[6,1] = widgets.Text(value='', description=dn, placeholder=ph, disabled=True)
        bt.cosmol2 = h, q, clambda, tu, ttg, zf, omega, omega_lambda

    out = widgets.interactive_output(f,{'o0':o0,'o1':o1,'o2':o2,'o3':o3,'o4':o4,'o5':o5,'o6':o6})

   # Create button widget to control flow
    button = widgets.Button(description='run zmag')
    button.style.button_color = 'lightblue'
    output = widgets.Output()
    grid[6,4] = button
    def f2(b):
        with output:
            if r1> 0:
                n = cf(bt.om)
            fplot(2,'3',1)
    display(grid, output)
    button.on_click(f2)

def namewidget(Model,Z,IMF,Mup,Library):
    # Builds file name according to dropdown menu selection
    # name = 'bc2020_z001_chab_lr_BaSeL_ssp.fits'
    # name = 'bc2020_z008_chab_MU300_hr_stelib_ssp.fits'

    filew = ''
    # Return if None selected
    if Model == 'None':
        return filew

    # Check for unavailable models
    if Model != 'CB19':
        Mup = '100'
        if 'Vaz' in IMF:
            IMF = 'Chabrier'

    # Check for undefined Z
    if 'None' in str(Z):
        Z = '0.017'

    if Model == 'BC03':
        filew = 'bc2003'
    elif Model == 'CB03':
        filew = 'cb2003'
    elif Model == 'CB07':
        filew = 'cb2007'
    elif Model == 'BC19':
        filew = 'bc2019'
    elif Model == 'CB19':
        filew = 'cb2019'
    elif Model == 'BC22':
        filew = 'bc2022'
    elif Model == 'CB22':
        filew = 'cb2022'

    if Z == '0':
        filew = filew + '_z0000_'
    elif Z == '0.0001':
        filew = filew + '_z0001_'
    elif Z == '0.0002':
        filew = filew + '_z0002_'
    elif Z == '0.0004':
        filew = filew + '_z0004_'
    elif Z == '0.0005':
        filew = filew + '_z0005_'
    elif Z == '0.001':
        filew = filew + '_z001_'
    elif Z == '0.002':
        filew = filew + '_z002_'
    elif Z == '0.004':
        filew = filew + '_z004_'
    elif Z == '0.006':
        filew = filew + '_z006_'
    elif Z == '0.008':
        filew = filew + '_z008_'
    elif Z == '0.010':
        filew = filew + '_z010_'
    elif Z == '0.014':
        filew = filew + '_z014_'
    elif Z == '0.017':
        filew = filew + '_z017_'
    elif Z == '0.019':
        filew = filew + '_z019_'
    elif Z == '0.02':
        filew = filew + '_z020_'
    elif Z == '0.03':
        filew = filew + '_z030_'
    elif Z == '0.04':
        filew = filew + '_z040_'
    elif Z == '0.05':
        filew = filew + '_z050_'
    elif Z == '0.06':
        filew = filew + '_z060_'
    elif Z == '0.100':
        filew = filew + '_z100_'

    if IMF == 'Chabrier':
        filew = filew + 'chab_'
    elif IMF == 'Kroupa':
        filew = filew + 'kroup_'
    elif IMF == 'Salpeter':
        filew = filew + 'salp_'
    elif IMF == 'Vaz_x0.3*':
        filew = filew + 'v0p30_'
    elif IMF == 'Vaz_x0.8*':
        filew = filew + 'v0p80_'
    elif IMF == 'Vaz_x1.0*':
        filew = filew + 'v1p00_'
    elif IMF == 'Vaz_x1.3*':
        filew = filew + 'v1p30_'
    elif IMF == 'Vaz_x1.5*':
        filew = filew + 'v1p50_'
    elif IMF == 'Vaz_x1.8*':
        filew = filew + 'v1p80_'
    elif IMF == 'Vaz_x2.0*':
        filew = filew + 'v2p00_'
    elif IMF == 'Vaz_x2.3*':
        filew = filew + 'v2p30_'
    elif IMF == 'Vaz_x2.8*':
        filew = filew + 'v2p80_'
    elif IMF == 'Vaz_x3.3*':
        filew = filew + 'v3p30_'

    if Mup == '10*':
        filew = filew + 'MU010_'
    elif Mup == '300*':
        filew = filew + 'MU300_'
    elif Mup == '600*':
        filew = filew + 'MU600_'

    if Library == 'Miles':
        filew = filew + 'hr_xmilesi_ssp'
    elif Library == 'Miles+':
        filew = filew + 'er_xmilesi_ssp'
    elif Library == 'IndoUS':
        filew = filew + 'hr_xindous_ssp'
    elif Library == 'Stelib':
        filew = filew + 'hr_stelib_ssp'
    elif Library == 'BaSeL':
        filew = filew + 'lr_BaSeL_ssp'

    return filew

def pl():
    # Executes function selected in widget1
    ff(nn)
    return

def pyGALAXEV():
    # For Jupyter Notebook

    # Define buffers
    global om,rm,nn,ff,save,jupy,cry
    bt.om    = ['']*12
    rm       = ['']*12
    om[4]    = 'None'
    om[5]    = 'None'
    bt.om[4] = om[4]
    bt.om[5] = om[5]
    nn = 0
    save    = False
    bt.jupy = True
    bt.ss   = False    # Don't store sed plots
    bt.sc   = False    # Don't store other plots

    # Define values to appear in file name Dropdown menus
    cb = {'None':['0.017','0','0.0001','0.0002','0.0005','0.001','0.002','0.004','0.006','0.008','0.010','0.014','0.017','0.02','0.03','0.04','0.06'],
          'BC19':['0.017','0','0.0001','0.0002','0.0005','0.001','0.002','0.004','0.006','0.008','0.010','0.014','0.017','0.02','0.03','0.04','0.06'],
          'CB19':['0.017','0','0.0001','0.0002','0.0005','0.001','0.002','0.004','0.006','0.008','0.010','0.014','0.017','0.02','0.03','0.04','0.06'],
         #'BC22':['0.017','0.0001','0.0005','0.001','0.002','0.004','0.006','0.008','0.010','0.014','0.017','0.02','0.03','0.04'],
         #'CB22':['0.017','0.0001','0.0005','0.001','0.002','0.004','0.006','0.008','0.010','0.014','0.017','0.02','0.03','0.04'],
          'BC03':['0.02','0.0001','0.0004','0.004','0.008','0.02','0.05'],
          'CB03':['0.019','0.0004','0.001','0.004','0.008','0.019','0.03'],
          'CB07':['0.02','0.0001','0.0004','0.004','0.008','0.02','0.05']}
    imf = ['Chabrier', 'Kroupa', 'Salpeter', 'Vaz_x0.3*', 'Vaz_x0.8*', 'Vaz_x1.0*', 'Vaz_x1.3*', 'Vaz_x1.5*', 'Vaz_x1.8*', 'Vaz_x2.0*', 'Vaz_x2.3*', 'Vaz_x2.8*', 'Vaz_x3.3*']
    mup = ['10*','100','300*','600*']
    atl = ['Miles','Miles+','IndoUS','Stelib','BaSeL']
    zW  = cb['CB19']
    dup = r'\(M_{UP}\)'

    # Define values to appear in SED plotting Dropdown menu
    vplt = ['SEDs', 'Magnitudes', 'Colors', 'Line indices', 'Physical Properties']
    # Define values to appear in task Dropdown menu
    vcod = ['None', 'csp_galaxev', 'add_bursts', 'cm_evolution', 'zmag', 'rf_phot', 'of_phot', 'bcfits2txt', 'galaxevpl']
    # Define default SEDs for model file names
    vsed = ['cb2019_z017_chab_hr_xmilesi_ssp', 'cb2019_z017_chab_hr_xmilesi_ssp']

    # Define widget grid
    grid = widgets.GridspecLayout(8,11)

    # Define copyright + revision date
    cry = '<b>pyGALAXEV</b>, Python Version (C) 2020-2023 - G. Bruzual and S. Charlot - All Rights Reserved (Last revision: ' + bs.fdate(os.environ.get('glxpyl') + '/.tar.timestamp',1) + ')'
    grid[0,0:10] = widgets.HTML(value=cry, placeholder='', description='')

    # Minimal instructions
    ins = '<span style="color:blue">Use the dropdown menus to build the file name corresponding to the desired model(s). Values marked (*) are available only for the <b>CB19</b> models. The command <b>model ID</b> provides details on model ingredients.</span>'
    grid[1,0:10] = widgets.HTML(value=ins, placeholder='', description='')

    # Dropdown menu's to build model file name 1
    grid[2,1] = widgets.Dropdown(options = cb.keys(), description='Model', value='CB19',layout={'width': 'max-content'})
    grid[2,2] = widgets.Dropdown(description='Z',layout={'width': 'max-content'})
    grid[2,3] = widgets.Dropdown(options=imf, description='IMF', value=imf[0], rows=len(imf), interactive=True, layout={'width': 'max-content'})      # Create dropdown widget to select IMF
    grid[2,4] = widgets.Dropdown(options=mup, description=dup, value=mup[1], rows=len(mup), interactive=True, layout={'width': 'max-content'})      # Create dropdown widget to select Mup
    grid[2,5] = widgets.Dropdown(options=atl, description='Library', value=atl[0], rows=len(atl), interactive=True, layout={'width': 'max-content'}) # Create dropdown widget to select Spectral Library
    mW1  = grid[2,1]
    zW1  = grid[2,2]
    imf1 = grid[2,3]
    mup1 = grid[2,4]
    lib1 = grid[2,5]
    def fn(Model, Z, IMF, Mup, Library):
        zW1.options = cb[Model] # Here is the trick, i.e. update zW.options based on model, namely modelW.value.
        fw = namewidget(Model,zW1.value,IMF,Mup,Library)
        bt.fw = fw
       #Textbox widget to select file
        grid[2,0] = widgets.Text(value=fw, description='File 1', disabled=False)
    out = widgets.interactive_output(fn,{'Model':mW1, 'Z':zW1, 'IMF':imf1, 'Mup':mup1, 'Library':lib1})

    # Dropdown menu's to build model file name 2
    grid[3,1] = widgets.Dropdown(options = cb.keys(), description='Model', value='CB19',layout={'width': 'max-content'})
    grid[3,2] = widgets.Dropdown(description='Z',layout={'width': 'max-content'})
    grid[3,3] = widgets.Dropdown(options=imf, description='IMF', value=imf[0], rows=len(imf), interactive=True, layout={'width': 'max-content'})      # Create dropdown widget to select IMF
    grid[3,4] = widgets.Dropdown(options=mup, description=dup, value=mup[1], rows=len(mup), interactive=True, layout={'width': 'max-content'})      # Create dropdown widget to select Mup
    grid[3,5] = widgets.Dropdown(options=atl, description='Library', value=atl[0], rows=len(atl), interactive=True, layout={'width': 'max-content'}) # Create dropdown widget to select Spectral Library
    mW2  = grid[3,1]
    zW2  = grid[3,2]
    imf2 = grid[3,3]
    mup2 = grid[3,4]
    lib2 = grid[3,5]
    def gn(Model, Z, IMF, Mup, Library):
        zW2.options = cb[Model] # Here is the trick, i.e. update zW.options based on model, namely modelW.value.
        gw = namewidget(Model,zW2.value,IMF,Mup,Library)
        bt.gw = gw
        # Textbox widget to select file 1
        grid[3,0] = widgets.Text(value=gw, description='File 2', disabled=False)
    out = widgets.interactive_output(gn,{'Model':mW2, 'Z':zW2, 'IMF':imf2, 'Mup':mup2, 'Library':lib2})

    # Minimal instructions
    ins = 'select the quantity to <b>Plot</b> or the <b>Task</b> to run and then click the <b>green</b> button. Wait for the file(s) to read. Use the <b>yellow</b> buttons for model information and the <b>cyan</b> buttons to create the indicated file.'
    grid[5,0:10] = widgets.HTML(value=ins, placeholder='', description='')
    #ins = 'Type or paste the model <b>file name(s)</b> in the space(s) provided, select the desired <b>task</b> and click the <b>green</b> or <b>yellow</b> buttons. Wait for the file(s) to read.'
    #grid[2,0:10] = widgets.HTML(value=ins, placeholder='', description='')
    #ins = '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<font size="1">Model files are searched for first in the current directory <span style="color:blue">./</span>, then in the distributed SSP model directory <span style="color:blue">glxssp</span>, and then in the output directory <span style="color:blue">glxout</span>. Output from codes you run is written to <span style="color:blue">glxout</span> and temporary files to <span style="color:blue">glxtmp</span>.</font>'
    #grid[3,0:10] = widgets.HTML(value=ins, placeholder='', description='')

    # Create dropdown widget to select plot
    grid[6,0] = widgets.Dropdown(options=vplt, description='Plot', value=vplt[0], rows=len(vplt), interactive=True)

    # Create button widget to select list models
    grid[6,2] = widgets.Button(description='Available Models')
    button = grid[6,2]
    button.style.button_color = 'yellow'
    output = widgets.Output()
    def f4(b):
        with output:
            bs.ldir ('l')
    button.on_click(f4)

    # Create button widget to select filter list
    grid[6,3] = widgets.Button(description='Filter List')
    button = grid[6,3]
    button.style.button_color = 'yellow'
    output = widgets.Output()
    def f5(b):
        with output:
            fl.filterid(0)
    button.on_click(f5)

    # Create button widget to select fits2ised
    grid[6,4] = widgets.Button(description='fits2ised')
    button = grid[6,4]
    button.style.button_color = 'cyan'
    output = widgets.Output()
    def f3(b):
        with output:
            om[3] = 'fits2ised'
            n, o, q = menu(1,0)
            fplot(n,o,q)
            om[3] = ''
    button.on_click(f3)

    # Create dropdown widget to select code to run
    grid[7,0] = widgets.Dropdown(options=vcod, description='Task', value=vcod[0], rows=len(vcod), interactive=True)

    # Create button widget to select model ID
    grid[7,2] = widgets.Button(description='Model ID')
    button = grid[7,2]
    button.style.button_color = 'yellow'
    output = widgets.Output()
    def f1(b):
        with output:
            o0 = grid[2,0]
            def fx(o0):
                global om
                om[0] = o0.value
                return om
            om = fx(o0)
            oo = [om[0],'q']
            bc.bcf2t(oo)
    button.on_click(f1)

    # Create button widget to select Release Notes
    grid[7,3] = widgets.Button(description='Release Notes')
    button = grid[7,3]
    button.style.button_color = 'yellow'
    output = widgets.Output()
    def f6(b):
        with output:
            bs.myhelp()
            bs.ldir('?')
    button.on_click(f6)

    # Create button widget to select bcfits2txt
    grid[7,4] = widgets.Button(description='bcfits2txt')
    button = grid[7,4]
    button.style.button_color = 'cyan'
    output = widgets.Output()
    def f9(b):
        global file1
        with output:
            o0 = grid[2,0]
            def fx(o0):
                global om
                om[0] = o0.value
                return om
            om = fx(o0)
            oo = [om[0],'q']
            bc.bcf2t(oo)
            file1 = om[0]
            widgetbcf2t()
    button.on_click(f9)

    # Create button widget to control flow
    button = widgets.Button(description='Plot / Run task')
    button.style.button_color = 'lime'
    output = widgets.Output()
    grid[7,1] = button
    def f7(b):
        global om
        bt.ss = False
        bt.sc = False
        with output:
            o0 = grid[2,0]
            o1 = grid[3,0]
            o2 = grid[6,0]
            o3 = grid[7,0]
            def f8(o0,o1,o2,o3):
                global om
                om[0] = o0.value
                om[1] = o1.value
                om[2] = o2.value
                om[3] = o3.value
                if om[1] == om[0]:
                    om[1] = ''
                return om
            om = f8(o0,o1,o2,o3)
            m = rf(1,om)-1
            n, o, q = menu(1,0)
            if o=='0':
                widgetcsp()
            elif o=='1':
                widgetadd()
            elif o=='2':
                widgetcmev()
            elif o=='3':
                widgetzmag()
            elif o=='4':
                widgetgpl()
            elif o=='5':
                widgetrfphot()
            elif o=='6':
                widgetofphot()
            elif o=='7':
                widgetbcf2t()
            else:
                gg(o,m)
                fplot(n,o,q)
    display(grid, output)
    button.on_click(f7)

def pinit(a):
    global om,jupy

    # Init parameters
    om = ['']*12
    om[5]  = 'None'
    om[6]  = 'None'
    om[10] = 'False'
    bt.jupy = False

    # Print header
    bs.head()
    bs.ldir('l')

    # get model filenames from command line
    n = len(a)
    if n==1:
        # Name of 1 file has been entered in command line
        om[0] = a[0]
        om[1] = ''
        n = rf(0,om)
    elif n==2:
        # Name of 2 files have been entered in command line
        om[0] = a[0]
        om[1] = a[1]
        n = rf(0,om)
    else:
        n=0
    return n

def cf(om):
    # Read model files
    global w1,f1,t1,h1,m1,d1,p1,a1,c1,e1,k1,qf1,file1

    # Check file names entered by user
    if om[0] == '':
        print('No file name has been entered')
    n=0
    file1=''
    if om[0] != '':
        n=n+1
        file1 = str(om[0])
        if file1.find('z') == 0:
            file1 = bs.zrep(bt.fw,file1)
            om[0] = file1
        file1 = bs.fcheck(file1)

    # Read fits table
    bt.f1 = bc.bcfits(file1)
    w1,f1,t1,e1,c1,m1,p1,d1,v1,a1,h1 = bt.f1
    k1 = 0
    qf1 = True
    hd.multhead(0,file1,t1,w1,e1)
    return n

def rf(k,om):
    # Read model files
    global w1,f1,t1,h1,m1,d1,p1,a1,c1,e1,k1,qf1,file1
    global w2,f2,t2,h2,m2,d2,p2,a2,c2,e2,k2,qf2,file2

    # Check file names entered by user
    if om[0] == '' and om[1] == '':
        print('No file name has been entered')
    n=1
    file1=''
    file2=''
    if om[0] != '':
        n=n+1
        file1 = str(om[0])
        if file1.find('z') == 0:
            file1 = bs.zrep(bt.fw,file1)
            om[0] = file1
        file1 = bs.fcheck(file1)
    if om[1] != '':
        n=n+1
        file2 = str(om[1])
        if file2.find('z') == 0:
            file2 = bs.zrep(bt.fw,file2)
            om[1] = file2
        file2 = bs.fcheck(file2)
        if len(file1) <= 0:
            file1 = file2
            file2 = ''
            n = 2

    # Choose correct file format for file 1
    if om[4] != file1:
        if '.fits' in file1:
            # Read fits table
            bt.f1 = bc.bcfits(file1)
            w1,f1,t1,e1,c1,m1,p1,d1,v1,a1,h1 = bt.f1
            k1 = 0
            qf1 = True
           #hd.multhead(0,file1,t1,w1,e1)
        elif '.ised' in file1:
            # Read .ised file
            w1,f1,t1,h1 = bc.read_bcised(file1)
            k1 = 1
            qf1 = False
        elif '.ased' in file1:
            # Read .ased file or any ascii table (galaxevpl output)
            w1,f1,t1,h1 = bc.read_bcased(file1)
            k1 = 1
            qf1 = True
        om[4] = file1

    # Choose correct file format for file 2
    if om[5] != file2:
        if '.fits' in file2:
            # Read fits table
            bt.f2 = bc.bcfits(file2)
            w2,f2,t2,e2,c2,m2,p2,d2,v2,a2,h2 = bt.f2
            k2 = 0
            qf2 = True
           #hd.multhead(0,file2,t2,w2,e2)
        elif '.ised' in file2:
            # Read .ised file
            w2,f2,t2,h2 = bc.read_bcised(file2)
            k2 = 1
            qf2 = False
        elif '.ased' in file2:
            # Read .ased file or any ascii table (galaxevpl output)
            w2,f2,t2,h2 = bc.read_bcased(file2)
            k2 = 1
            qf2 = True
        om[5] = file2
    return n

def bands():
    # Print photometric bands available in fits file

    # '         Mbol', '    U_Johnson', '   B2_Johnson', '   B3_Johnson', '    V_Johnson', '    R_Johnson', '    I_Johnson', '    J_Johnson', '    K_Johnson', '    L_Johnson',
    # '    R_Cousins', '    I_Cousins', '    J_Palomar', '    H_Palomar', '    K_Palomar', '      J_2Mass', '      H_2Mass', '     Ks_2Mass', ' Kprime_Cowie', '    I3p6_IRAC',
    # '    I4p5_IRAC', '    I5p7_IRAC', '    I7p9_IRAC', '     I12_IRAS', '     I25_IRAS', '     I60_IRAS', '    I100_IRAS', '     M24_MIPS', '     M70_MIPS', '    M160_MIPS',
    # '       u_SDSS_AB', '       g_SDSS_AB', '       r_SDSS_AB', '       i_SDSS_AB', '       z_SDSS_AB', '   u1_CFHT_MC_AB', '   u3_CFHT_MC_AB', '   g1_CFHT_MC_AB', '   g3_CFHT_MC_AB', '   r1_CFHT_MC_AB',
    # '   r3_CFHT_MC_AB', '   i2_CFHT_MC_AB', '   i3_CFHT_MC_AB', '   z1_CFHT_MC_AB', '   z3_CFHT_MC_AB', '      H_2MASS_AB', '   Ks_CFHT_WC_AB', '    FUV_GALEX_AB', '    NUV_GALEX_AB',
    # '   WFC3_F225W_AB', '   WFC3_F336W_AB', '  WFC3_FR388N_AB', '   WFC3_F438W_AB', '   WFC3_F555W_AB', '   WFC3_F814W_AB', '   WFC3_F110W_AB', '   WFC3_F125W_AB', '   WFC3_F160W_AB',
    # '  UVIS1_f225w_AB', '  UVIS1_f275w_AB', '  UVIS1_f336w_AB', '  UVIS1_f438w_AB', '  UVIS1_f555w_AB', '  UVIS1_f547m_AB', '  UVIS1_f606w_AB', '  UVIS1_f625w_AB', '  UVIS1_f656n_AB',
    # '  UVIS1_f657n_AB', '  UVIS1_f658n_AB', '  UVIS1_f814w_AB', ' ACSWFC_F220w_AB', ' ACSWFC_F250w_AB', ' ACSWFC_F330w_AB', ' ACSWFC_F410w_AB', ' ACSWFC_F435w_AB', ' ACSWFC_F475w_AB',
    # ' ACSWFC_F555w_AB', ' ACSWFC_F606w_AB', ' ACSWFC_F625w_AB', ' ACSWFC_F775w_AB', ' ACSWFC_F814w_AB',
    b = ['Mbol','U_Johnson','B2_Johnson','B3_Johnson','V_Johnson','R_Johnson','I_Johnson','J_Johnson','K_Johnson','L_Johnson','R_Cousins','I_Cousins','J_Palomar','H_Palomar','K_Palomar',
         'J_2Mass','H_2Mass','Ks_2Mass','Kprime_Cowie','I3p6_IRAC','I4p5_IRAC','I5p7_IRAC','I7p9_IRAC','I12_IRAS','I25_IRAS','I60_IRAS','I100_IRAS','M24_MIPS','M70_MIPS','M160_MIPS',
         'u_SDSS_AB','g_SDSS_AB','r_SDSS_AB','i_SDSS_AB','z_SDSS_AB','u1_CFHT_MC_AB','u3_CFHT_MC_AB','g1_CFHT_MC_AB','g3_CFHT_MC_AB','r1_CFHT_MC_AB','r3_CFHT_MC_AB','i2_CFHT_MC_AB',
         'i3_CFHT_MC_AB','z1_CFHT_MC_AB','z3_CFHT_MC_AB','H_2MASS_AB','Ks_CFHT_WC_AB','FUV_GALEX_AB','NUV_GALEX_AB','WFC3_F225W_AB','WFC3_F336W_AB','WFC3_FR388N_AB','WFC3_F438W_AB',
         'WFC3_F555W_AB','WFC3_F814W_AB','WFC3_F110W_AB','WFC3_F125W_AB','WFC3_F160W_AB','UVIS1_f225w_AB','UVIS1_f275w_AB','UVIS1_f336w_AB','UVIS1_f438w_AB','UVIS1_f555w_AB','UVIS1_f547m_AB',
         'UVIS1_f606w_AB','UVIS1_f625w_AB','UVIS1_f656n_AB','UVIS1_f657n_AB','UVIS1_f658n_AB','UVIS1_f814w_AB','ACSWFC_F220w_AB','ACSWFC_F250w_AB','ACSWFC_F330w_AB','ACSWFC_F410w_AB',
         'ACSWFC_F435w_AB','ACSWFC_F475w_AB','ACSWFC_F555w_AB','ACSWFC_F606w_AB','ACSWFC_F625w_AB','ACSWFC_F775w_AB','ACSWFC_F814w_AB']
    print()
    print('Photometric bands available in this fits file:')
    print()
    print('  N        Band      N       Band      N           Band      N            Band      N             Band')
    print('  -------------     -------------     -----------------     ------------------     -------------------')
    print('  0:       Mbol     15:   J_2Mass     30:     u_SDSS_AB     49:  WFC3_F225W_AB     70: ACSWFC_F220w_AB')
    print('                    16:   H_2Mass     31:     g_SDSS_AB     50:  WFC3_F336W_AB     71: ACSWFC_F250w_AB')
    print('  1:  U_Johnson     17:  Ks_2Mass     32:     r_SDSS_AB     51: WFC3_FR388N_AB     72: ACSWFC_F330w_AB')
    print('  2: B2_Johnson                       33:     i_SDSS_AB     52:  WFC3_F438W_AB     73: ACSWFC_F410w_AB')
    print("  3: B3_Johnson     18:  K'_Cowie     34:     z_SDSS_AB     53:  WFC3_F555W_AB     74: ACSWFC_F435w_AB")
    print('  4:  V_Johnson                                             54:  WFC3_F814W_AB     75: ACSWFC_F475w_AB')
    print('  5:  R_Johnson     19: I3.6_IRAC     35: u1_CFHT_MC_AB     55:  WFC3_F110W_AB     76: ACSWFC_F555w_AB')
    print('  6:  I_Johnson     20: I4.5_IRAC     36: u3_CFHT_MC_AB     56:  WFC3_F125W_AB     77: ACSWFC_F606w_AB')
    print('  7:  J_Johnson     21: I5.7_IRAC     37: g1_CFHT_MC_AB     57:  WFC3_F160W_AB     78: ACSWFC_F625w_AB')
    print('  8:  K_Johnson     22: I7.9_IRAC     38: g3_CFHT_MC_AB                            79: ACSWFC_F775w_AB')
    print('  9:  L_Johnson                       39: r1_CFHT_MC_AB     58: UVIS1_f225w_AB     80: ACSWFC_F814w_AB')
    print('                    23:  I12_IRAS     40: r3_CFHT_MC_AB     59: UVIS1_f275w_AB')
    print(' 10:  R_Cousins     24:  I25_IRAS     41: i2_CFHT_MC_AB     60: UVIS1_f336w_AB')
    print(' 11:  I_Cousins     25:  I60_IRAS     42: i3_CFHT_MC_AB     61: UVIS1_f438w_AB')
    print('                    26: I100_IRAS     43: z1_CFHT_MC_AB     62: UVIS1_f555w_AB')
    print(' 12:  J_Palomar                       44: z3_CFHT_MC_AB     63: UVIS1_f547m_AB')
    print(' 13:  H_Palomar     27:  M24_MIPS     45:    H_2MASS_AB     64: UVIS1_f606w_AB')
    print(' 14:  K_Palomar     28:  M70_MIPS     46: Ks_CFHT_WC_AB     65: UVIS1_f625w_AB')
    print('                    29: M160_MIPS                           66: UVIS1_f656n_AB')
    print('                                      47:  FUV_GALEX_AB     67: UVIS1_f657n_AB')
    print('                                      48:  NUV_GALEX_AB     68: UVIS1_f658n_AB')
    print('                                                            69: UVIS1_f814w_AB')
    print(' -----------------------------------------------------------------------------------------------------')
    return b

def col(n):
    # Plot colors for option n

    def ucol_1(m):
        global bx
        # Plot user selected colors
        if not bt.jupy:
            # Ask for colors to plot, 1 model
            wg = input('Select up to 5 bands to plot up to 4 colors: N1 N2 N3 N4 N5 = ')
            wg = wg.replace(",", " ")
            wg = wg.split()
        else:
            wg = om[7]
        n = min(5,len(wg))
        f = [""]*n
        for i in range(n):
            f[i] = bx[int(wg[i])]
        k = n-1
        if k <= 0:
            p=0.
            l=''
            n=-1
            return n,p,l
        # k colors to plot
        p = [0.]*k
        l = [""]*k
        for i in range (k):
            p[i] = m[f[i]]-m[f[i+1]]
            l[i] = f[i] + ' - ' + f[i+1]
            l[i] = l[i].replace("prime", "'")
            l[i] = l[i].replace("p", ".")
        return k,p,l

    def ucol_2(m1,m2):
        global bx
        # Plot user selected colors
        if not bt.jupy:
            # Ask for colors to plot, 2 models
            wg = input('Select up to 5 bands to plot up to 4 colors: N1 N2 N3 N4 N5 = ')
            wg = wg.replace(",", " ")
            wg = wg.split()
        else:
            wg = om[7]
        n = min(5,len(wg))
        f = [""]*n
        for i in range(n):
            f[i] = bx[int(wg[i])]
        k = n-1
        if k <= 0:
            p1=0.
            p2=0.
            l=''
            n=-1
            return n,p1,p2,l
        # k colors to plot
        p1 = [0.]*k
        p2 = [0.]*k
        l  = [""]*k
        for i in range (k):
            p1[i] = m1[f[i]]-m1[f[i+1]]
            p2[i] = m2[f[i]]-m2[f[i+1]]
            l[i]  = f[i] + ' - ' + f[i+1]
            l[i]  = l[i].replace("prime", "'")
            l[i]  = l[i].replace("p", ".")
        return k,p1,p2,l

    def col_h(n):
        # Returns header and labels for color and magnitude plots
        if (n == 1):
            lh = ['U_Johnson','B2_Johnson',  'B3_Johnson','V_Johnson',  'V_Johnson','R_Johnson',
                'R_Johnson','I_Johnson',   'I_Johnson','J_Johnson',   'J_Johnson','K_Johnson',
                'K_Johnson','L_Johnson',   'V_Johnson','J_Johnson',   'V_Johnson','K_Johnson']
            ly = ['U - B', 'B - V', 'V - R', 'R - I', 'I - J', 'J - K', 'K - L', 'V - J', 'V - K']
        elif (n==2):
            lh = ['U_Johnson','B2_Johnson',  'B3_Johnson','V_Johnson',  'V_Johnson','R_Cousins',
                'R_Cousins','I_Cousins',   'I_Cousins','J_Johnson',   'J_Johnson','K_Johnson',
                'K_Johnson','L_Johnson',   'V_Johnson','J_Johnson',   'V_Johnson','K_Johnson']
            ly = ['U - B', 'B - V', 'V - R$_c$', 'R$_c$ - I$_c$', 'I$_c$ - J', 'J - K', 'K - L', 'V - J', 'V - K']
        elif (n==3):
            lh = ['u_SDSS_AB','g_SDSS_AB',   'g_SDSS_AB','r_SDSS_AB',   'r_SDSS_AB','i_SDSS_AB',   'i_SDSS_AB','z_SDSS_AB']
            ly = ['u - g', 'g - r', 'r - i', 'i - z']
        return lh,ly

    print()
    s=n
    if (n==0):
        n,f1,l = ucol_1(m1)
        if n > 0:
            nplt(n,a1,f1,l,file1,False,save)

    elif (n<=2):
        lh,ly = col_h(n)
        splt(3,3,a1,m1,lh,ly,file1,False,save)

    elif (n==3):
        lh,ly = col_h(n)
        splt(2,2,a1,m1,lh,ly,file1,False,save)

    elif (n==4):
        n,f1,f2,l = ucol_2(m1,m2)
        if n > 0:
            nplt2(n,a1,f1,a2,f2,l,file1,file2,False,save)

    elif (n<=6):
        lh,ly = col_h(n-4)
        splt2(3,3,a1,m1,a2,m2,lh,ly,file1,file2,False,save)

    elif (n==7):
        lh,ly = col_h(n-4)
        splt2(2,2,a1,m1,a2,m2,lh,ly,file1,file2,False,save)

    if n > 0 and bt.sc:
        # Rename saved file
        if s <= 4:
            filout = os.environ.get('glxout') + '/' + bs.lnam(file1).replace('fits','color') + str(s) + '.png'
        else:
            filout = os.environ.get('glxout') + '/' + bs.lnam(file1).replace('fits','') + bs.lnam(file2).replace('fits','color') + str(s-4) + '.png'
        filtmp = os.environ.get('glxtmp') + '/' + 'tmpplot.png'
        os.rename(filtmp,filout)
        print(bs.gp(' Plot stored in file: ' + filout))
    return n

def fcol(off):
    # Plot selected color sets
    global bx
    if not bt.jupy:
        bx = bands()
        while True:
            # Ask for colors to plot
            print()
            print('Select colors to plot (press <Enter> to exit loop):')
            print('        -1: UBVRIJKL - Johnson')
            print('        -2: UBVRIJKL - Johnson-Cousins')
            print('        -3:    ugriz - SDSS')
            print('        -4: user selected colors')
            op = input( 'choice = ')
            if len(op) <= 0 or op == 'q':
               break
            op=gop(op,off,4)
            col(op)
    else:
        # Show available bands
        cv = ['UBVRIJKL - Johnson', 'UBVRIJKL - Johnson-Cousins', 'ugriz - SDSS']
        bx = bands()
        d1 = 'Preselected'
        d2 = 'Bands'
        ph = 'Enter up to 5 bands'
        widget1(cv,d1,d2,ph,off,col)

def fidx(off):
    # Plot selected line indices

    def indices():
        # Print line strength indices available in fits file

        # 'logage'   , 'CN_1     ', 'CN_2     ', 'Ca4227   ', 'G4300    ', 'Fe4383   ', 'Ca4455   ', 'Fe4531   ', 'Fe4668   ', 'Hbeta    ', 'Fe5015   ', 'Mg_1     ', 'Mg_2     ', 'Mg-b     ', 'Fe5270   ',
        # 'Fe5335   ', 'Fe5406   ', 'Fe5709   ', 'Fe5782   ', 'Na-D     ', 'TiO_1    ', 'TiO_2    ', 'HdeltaA  ', 'HgammaA  ', 'HdeltaF  ', 'HgammaF  ', 'D4000    ', 'B4000_VN ', 'CaII8498 ',
        # 'CaII8542 ', 'CaII8662 ', 'MgI8807  ', 'H8_3889  ', 'H9_3835  ', 'H10_3798 ', 'BH-HK    ', 'BL1302   ', 'SiIV     ', 'BL1425   ', 'Fe1453   ', 'CIV1548a ', 'CIV1548c ', 'CIV1548e ',
        # 'BL1617   ', 'BL1664   ', 'BL1719   ', 'BL1853   ', 'FeII2402 ', 'BL2538   ', 'FeII2609 ', 'MgII     ', 'MgI      ', 'Mgwide   ', 'FeI      ', 'BL3096   ', 'CIVabs   ', 'HeIIems  '
        b = ['logage','CN_1','CN_2','Ca4227','G4300','Fe4383','Ca4455','Fe4531','Fe4668','Hbeta','Fe5015','Mg_1','Mg_2','Mg-b','Fe5270',
            'Fe5335','Fe5406','Fe5709','Fe5782','Na-D','TiO_1','TiO_2','HdeltaA','HgammaA','HdeltaF','HgammaF','D4000','B4000_VN','CaII8498',
            'CaII8542','CaII8662','MgI8807','H8-3889','H9-3835','H10-3798','BH-HK','BL1302','SiIV','BL1425','Fe1453','CIV1548a','CIV1548c','CIV1548e',
            'BL1617','BL1664','BL1719','BL1853','FeII2402','BL2538','FeII2609','MgII','MgI','Mgwide','FeI','BL3096','CIVabs','HeIIems']
        u = ['yr','(mag)','(mag)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(mag)','(mag)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(mag)','(mag)','(A)','(A)','(A)','(A)','','','(A)','(A)','(A)',
            '(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)']
        print()
        print('Line strength indices available in this fits file:')
        print()
        print('  <----- Lick Indices ----->      <-- W&O -->      <--- DTT -->      <- Fanelli et al. UV Indices ->')
        print('  N  Index         N  Index        N  Index         N  Index          N  Index           N  Index')
        print('  --------------------------      -----------      ------------      -------------------------------')
        print('  1: CN_1         11: Mg_1        22: HdeltaA      28: CaII8498      36: BL1302         47: FeII2402')
        print('  2: CN_2         12: Mg_2        23: HgammaA      29: CaII8542      37: SiIV           48: BL2538')
        print('  3: Ca4227       13: Mg-b        24: HdeltaF      30: CaII8662      38: BL1425         49: FeII2609')
        print('  4: G4300        14: Fe5270      25: HgammaF      31: MgI8807       39: Fe1453         50: MgII')
        print('  5: Fe4383       15: Fe5335                                         40: CIV1548a       51: MgI')
        print('  6: Ca4455       16: Fe5406      <- 4000A ->      <Marcillac>       41: CIV1548c       52: Mgwide')
        print('  7: Fe4531       17: Fe5709      ----------       -----------       42: CIV1548e       53: FeI')
        print('  8: Fe4668       18: Fe5782      26: B4000_GC     32: H8-3889       43: BL1617         54: BL3096')
        print('  9: Hbeta        19: Na-D        27: B4000_VN     33: H9-3835       44: BL1664         55: CIVabs')
        print(' 10: Fe5015       20: TiO_1                        34: H10-3798      45: BL1719         56: HeIIems')
        print('                  21: TiO_2                        35: BH-HK         46: BL1853')
        print(' ---------------------------------------------------------------------------------------------------')
        return b

    global bx
    if not bt.jupy:
        bx = indices()
        while True:
            # Ask for indices to plot
            print()
            print('Select indices to plot (press <Enter> to exit loop):')
            print('        -1: Lick Indices  1:15')
            print('        -2: Lick Indices 16:21 + WO + Marcillac et al. indices')
            print('        -3: DTT + HK + B4000vn indices')
            print('        -4: Fanelli et al. UV indices')
            print('        -5: Fanelli et al. UV indices (cont.)')
            print('        -6: user selected indices')
            print('         N: to plot index N')
            op = input( 'choice = ')
            if len(op) <= 0 or op == 'q':
                break
            op=gop(op,off,6)
            idx(op)
    else:
        # Show available indices
        cv = ['Lick Indices  1:15', 'Lick Indices 16:21 + WO + Marcillac et al. indices', 'DTT + HK + B4000vn indices',
              'Fanelli et al. UV indices', 'Fanelli et al. UV indices (cont.)']
        bx = indices()
        d1 = 'Preselected'
        d2 = 'Indices'
        ph = 'Enter up to 4 indices'
        widget1(cv,d1,d2,ph,off,idx)

def fmag(off):
    # Plot selected photometric magnitudes
    global bx
    if not bt.jupy:
        bx = bands()
        while True:
            # Ask for magnitudes to plot
            print()
            print('Select magnitudes to plot (press <Enter> to exit loop):')
            print('        -1: UV to FIR')
            print('        -2: NIR to FIR')
            print('        -3: user selected bands')
            print('         N: to plot band N')
            op = input( 'choice = ')
            if len(op) <= 0 or op == 'q':
                break
            op=gop(op,off,3)
            mag(op)
    else:
        # Show available bands
        cv = ['UV to FIR', 'NIR to FIR']
        bx = bands()
        d1 = 'Preselected'
        d2 = 'Bands'
        ph = 'Enter up to 4 bands'
        widget1(cv,d1,d2,ph,off,mag)

def fphy(off):
    # Plots selected physical properties

    def properties():
        # Print physical properties available in fits file
        b = ['logage', 'logNLy', 'logNHeI', 'logNHeII', 'logNHeII-logNLy', 'B912', 'B4000VN', 'B4000SDSS', 'B4000', 'SNR', 'PISNR', 'RegSNR', 'TypeIaSNR', 'FailedSNR',
            'PNBR', 'NBH', 'NNS', 'NWD', 'Lx', 'Mstars', 'Mremnants', 'Mretgas', 'Mgalaxy', 'SFR', 'Mtot', 'Mtot_Lb', 'Mtot_Lv', 'Mtot_Lk', 'Mliv_Lb', 'Mliv_Lv', 'Mliv_Lk',
            'EvolutionaryFlux', 'SpecificEvolFlux', 'TurnoffMass', 'BPMS_BMS', 'TotalMassLossRate', 'DustProductionRate']
        u = ['(yr)','(photons/Mo)','(photons/Mo)','(photons/Mo)','','','','','', 'SN/yr', 'SN/yr', 'SN/yr', 'SN/yr', 'SN/yr', 'PN/yr', 'number/Mo', 'number/Mo', 'number/Mo', '(Lo)',
            '(Mo)','(Mo)', '(Mo)','(Mo)','(Mo/yr)','(Mo)','(Mo/Lo)','(Mo/Lo)','(Mo/Lo)','(Mo/Lo)','(Mo/Lo)','(Mo/Lo)','(Nstars/yr)','(Nstars/yr/Lo)','(Mo)','','(M/Mo/yr)','(M/Mo/yr)']
        print()
        print('Physical properties available in this fits file:')
        print()
        print('  N  Property                             N  Property                              N  Property')
        print('  ---------------------------------      -----------------------------------      ------------------------------')
        print('  1: HI ionizing photons                 14: Planetary nebula birth rate          24: Mtot = Mliv + Mrem')
        print('  2: HeI ionizing photons                15: Number of black holes                25: Mtot/L(B)')
        print('  3: HeII ionizing photons               16: Number of neutron stars              26: Mtot/L(V)')
        print('  4: Ratio HeII/HI ionizing photons      17: Number of white dwarfs               27: Mtot/L(K)')
        print('                                         18: X-ray binary luminosity              28: Mliv/L(B)')
        print('  5:  912 A break amplitude                                                       29: Mliv/L(V)')
        print('  6: 4000 A break amplitude (VN)         19: Mliv = Mass in living stars          30: Mliv/L(K)')
        print('  7: 4000 A break amplitude (SDSS)       20: Mrem = Mass in stellar remnants')
        print('  8: 4000 A break amplitude              21: Mgas = Mass in processed gas         31: Evolutionary flux')
        print('                                         22: Mgal = Mass of galaxy                32: Specific evolutionary flux')
        print('  9: Total supernova rate                         = Mliv+ Mrem + Mgas             33: Turnoff mass')
        print(' 10: Pair Instability supernova rate                                              34: Bol flux PMS / Bol flux MS')
        print(' 11: Regular supernova rate              23: Star formation rate                  35: Total mass loss rate')
        print(' 12: Type Ia supernova rate                                                       36: Dust production rate')
        print(' 13: Failed supernova rate')
        print(' ---------------------------------------------------------------------------------------------------------------')
        return b

    global bx
    if not bt.jupy:
        bx= properties()
        while True:
            # Ask for indices to plot
            print()
            print('Select quantities to plot (press <Enter> to exit loop):')
            print('        -1: Mass related properties')
            print('        -2: Ionizing photons and remnant production rate')
            print('        -3: user selected property')
            print('         N: to plot property N')
            op = input( 'choice = ')
            if len(op) <= 0 or op == 'q':
                break
            op=gop(op,off,3)
            phy(op)
    else:
        # Show available properties
        cv = ['Mass related properties', 'Ionizing photons and remnant production rate']
        bx = properties()
        d1 = 'Properties'
        d1 = 'Preselected'
        d2 = 'Properties'
        ph = 'Enter up to 4 properties'
        widget1(cv,d1,d2,ph,off,phy)

def fplot1(n,o):

    global h, omega, omega_lambda, clambda, q, tu, ttg, zf

    while True:
        # Select plot option
        if o=='s' or o=='n':
            wn = normwav(o)
            fsed(0)
        elif o=='c':
            fcol(1)
        elif o=='m':
            fmag(1)
        elif o=='i':
            fidx(1)
        elif o=='p':
            fphy(1)
        elif o == 'l':
            bs.ldir('l')
        elif o == 'h' or o=='?':
            bs.myhelp()
            bs.ldir('?')
        elif o=='f':
            n=0
        elif o=='0':
            if not bt.jupy:
                csp.csp_galaxev(file1)
            else:
                inp = wcsp(file1)
                csp.csp_galaxev(inp)
        elif o=='1':
            if not bt.jupy:
                add.pyadd_bursts(file1)
            else:
                inp = wadd()
                add.pyadd_bursts(inp)
        elif o=='2':
            if not bt.jupy:
                evl.cmev(file1,bt.jupy)
            else:
                inp = wcmev()
                evl.cmev(inp,bt.jupy)
        elif o=='3':
            if not bt.jupy:
                evl.zmag(file1)
            else:
                wzmag(file1)
        elif o=='4':
            if not bt.jupy:
                bc.glxpl(file1)
            else:
                bc.wglxpl(file1)
        elif o=='5':
            if not bt.jupy:
                evl.krfp(file1)
            evl.rf_phot()
        elif o=='6':
            if not bt.jupy:
                evl.kofp(file1)
            evl.of_phot()
        elif o=='7':
            if not bt.jupy:
                bc.bcf2t([file1])
            else:
                bc.bcf2t(bt.om)
        elif o=='10':
            hd.idheader(e1,len(w1))
        elif o=='11':
            fl.filterid(0)
        else:
            print('Unknown option: ' + o)
        return n

def fplot2(n,o):

    while True:
       # Select plot options for 2 files
        if o=='s' or o=='n':
            wn = normwav(o)
            fsed(4)
        elif o=='c':
            fcol(5)
        elif o=='m':
            fmag(4)
        elif o=='i':
            fidx(7)
        elif o=='p':
            fphy(4)
        elif o=='10':
            hd.idheader(e1,len(w1))
            hd.idheader(e2,len(w2))
        return n

def fsed(off):
    # Plots model SED at indicated age
    if not bt.jupy:
        # Ask for plot limits
        print()
        b = input(' Plot limits: xmin xmax log(ymin) log(ymax) [full range] = ')
        if len(b) > 0:
            b = b.replace(",", " ")
        om[8] = b.split()
        print()
        om[9] = input(' Normalize at wavelength (A) [none] = ')
        print()
        b = input(' Save fig(s) to png file(s) (y/n) [n] = ')
        if b == 'y' or b == 'Y':
            bt.ss = True
        else:
            bt.ss = False
        while True:
            # Ask for time steps to plot
            print()
            print(' Enter age of sed''s to plot, separated by space. Type ENTER to exit loop.')
            b = input('    (Hint: -t will select a set of spectra not exceeding age t) = ')
            if len(b) <= 0:
                break
            # Plot sed's
            b = b.replace(',', '|')
            b = b.replace(' ', '|')
            om[7] = b
            sed(off)
    else:
        # Show options
        widget2(off,sed)

def gop(op,off,m):
    # Read option
    global one
    op = -1*int(op)
    if op < 0:
        one = -op
        op  = m
    else:
        one = 0
    if op <m:
        op = op + off - 1
    else:
        if off < m:
           op = 0
        else:
           op = m
    return op

def idx(n):
    # Plot indices for option n

    def uidx_1(m):
        global bx, one
        # Plot user selected indices
        # Index units
        u = ['yr','(mag)','(mag)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(mag)','(mag)','(A)','(A)',
            '(A)','(A)','(A)','(A)','(A)','(mag)','(mag)','(A)','(A)','(A)','(A)','','','(A)','(A)','(A)','(A)',
            '(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)',
            '(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)']
        if not bt.jupy:
            if one==0:
                # Ask for indices to plot, 1 model
                wg = input('Select up to 4 indices to plot (<Enter> exits) N1 N2 N3 N4 = ')
                wg = wg.replace(",", " ")
                wg = wg.split()
                n = min(4,len(wg))
            else:
                wg = [str(one)]
                n = 1
        else:
            wg = om[7]
            n = min(4,len(wg))
        f = [""]*n
        for i in range(n):
            f[i] = bx[int(wg[i])]
        if n <= 0:
            p=0.
            l=''
            n=-1
            return n,p,l
        # k bands to plot
        p = [0.]*n
        l = [""]*n
        for i in range (n):
            l[i] = f[i] + ' ' + u[int(wg[i])]
            f[i] = f[i].replace("_", "")
            f[i] = f[i].replace("-", "")
            p[i] = m[f[i]]
        return n,p,l

    def uidx_2(m1,m2):
        global bx, one
        # Plot user selected indices
        # Index units
        u = ['yr','(mag)','(mag)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(mag)','(mag)','(A)','(A)',
            '(A)','(A)','(A)','(A)','(A)','(mag)','(mag)','(A)','(A)','(A)','(A)','','','(A)','(A)','(A)','(A)',
            '(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)',
            '(A)','(A)','(A)','(A)','(A)','(A)','(A)','(A)']
        if not bt.jupy:
            if one==0:
                # Ask for index to plot, 2 models
                wg = input('Select up to 4 indices to plot (<Enter> exits) N1 N2 N3 N4 = ')
                wg = wg.replace(",", " ")
                wg = wg.split()
                n = min(4,len(wg))
            else:
                wg = [str(one)]
                n = 1
        else:
            wg = om[7]
            n = min(4,len(wg))
        f = [""]*n
        for i in range(n):
            f[i] = bx[int(wg[i])]
        if n <= 0:
            p1=0.
            p2=0.
            l=''
            n=-1
            return n,p1,p2,l
        # k bands to plot
        p1 = [0.]*n
        p2 = [0.]*n
        l  = [""]*n
        for i in range (n):
            l[i]  = f[i] + ' ' + u[int(wg[i])]
            f[i]  = f[i].replace("_", "")
            f[i]  = f[i].replace("-", "")
            p1[i] = m1[f[i]]
            p2[i] = m2[f[i]]
        return n,p1,p2,l

    def idx_h(n):
        # Returns header and labels for index plots
        if (n == 1):
            # Lick indices 1:15
            l = ['CN1:m', 'CN2:m', 'Ca4227', 'G4300', 'Fe4383', 'Ca4455', 'Fe4531', 'Fe4668', 'Hbeta', 'Fe5015', 'Mg1:m', 'Mg2:m', 'Mgb', 'Fe5270', 'Fe5335']
        elif (n == 2):
            # Lick indices 16:21 + WO + Marcillac et al. indices
            l = ['Fe5406', 'Fe5709', 'Fe5782', 'NaD', 'TiO1:m', 'TiO2:m', 'HdeltaA', 'HgammaA', 'HdeltaF', 'HgammaF', '', '', 'H83889', 'H93835', 'H103798']
        elif (n==3):
            # DTT + HK + B4000 indices
            l = ['CaII8498', 'CaII8542', 'CaII8662', 'MgI8807', '', '', 'BHHK', '', '', 'B4000VN:n', '', '', '', '', '']
        elif (n==4):
            # Fanelli et al UV indices
            l = ['BL1302','SiIV','BL1425','Fe1453','CIV1548a','CIV1548c','CIV1548e','BL1617','BL1664','BL1719','BL1853','FeII2402','BL2538','FeII2609','MgII']
        elif (n==5):
            # Fanelli et al UV indices (cont.)
            l = ['MgI','Mgwide','FeI','BL3096','CIVabs','HeIIems','','','','','','','','','']
        lh,ly = yl(l)
        return lh,ly

    print()
    s=n
    if (n==0):
        n,f1,l = uidx_1(d1)
        if n > 0:
            nplt(n,a1,f1,l,file1,False,save)

    elif (n<=5):
        lh,ly = idx_h(n)
        splt(5,3,a1,d1,lh,ly,file1,False,save)

    elif (n==6):
        n,f1,f2,l = uidx_2(d1,d2)
        if n > 0:
            nplt2(n,a1,f1,a2,f2,l,file1,file2,False,save)

    elif (n<=11):
        lh,ly = idx_h(n-6)
        splt2(5,3,a1,d1,a2,d2,lh,ly,file1,file2,False,save)

    if n > 0 and bt.sc:
        # Rename saved file
        if s <=6:
            filout = os.environ.get('glxout') + '/' + bs.lnam(file1).replace('fits','lsindx') + str(s) + '.png'
        else:
            filout = os.environ.get('glxout') + '/' + bs.lnam(file1).replace('fits','') + bs.lnam(file2).replace('fits','lsindx') + str(s) + '.png'
        filtmp = os.environ.get('glxtmp') + '/' + 'tmpplot.png'
        os.rename(filtmp,filout)
        print(bs.gp(' Plot stored in file: ' + filout))
    return n

def mag(n):
    # Plot magnitudes for option n

    def umag_1(m):
        global bx, one
        # Plot user selected magnitudes
        if not bt.jupy:
            if one==0:
                # Ask for band to plot, 1 model
                wg = input('Select up to 4 bands to plot (<Enter> exits) N1 N2 N3 N4 = ')
                wg = wg.replace(",", " ")
                wg = wg.split()
                n = min(4,len(wg))
            else:
                wg = [str(one)]
                n = 1
        else:
            wg = om[7]
            n = min(4,len(wg))
        f = [""]*n
        for i in range(n):
            f[i] = bx[int(wg[i])]
        if n <= 0:
            p=0.
            l=''
            n=-1
            return n,p,l
        # k bands to plot
        p = [0.]*n
        l = [""]*n
        for i in range (n):
            p[i] = m[f[i]]
            l[i] = f[i]
            l[i] = l[i].replace("prime", "'")
            l[i] = l[i].replace("p", ".")
        return n,p,l

    def umag_2(m1,m2):
        global bx, one
        # Plot user selected magnitudes
        if not bt.jupy:
            if one==0:
                # Ask for band to plot, 1 model
                wg = input('Select up to 4 bands to plot (<Enter> exits) N1 N2 N3 N4 = ')
                wg = wg.replace(",", " ")
                wg = wg.split()
                n = min(4,len(wg))
            else:
                wg = [str(one)]
                n = 1
        else:
            wg = om[7]
            n = min(4,len(wg))
        f = [""]*n
        for i in range(n):
            f[i] = bx[int(wg[i])]
        if n <= 0:
            p1=0.
            p2=0.
            l=''
            n=-1
            return n,p1,p2,l
        # k bands to plot
        p1 = [0.]*n
        p2 = [0.]*n
        l  = [""]*n
        for i in range (n):
            p1[i] = m1[f[i]]
            p2[i] = m2[f[i]]
            l[i]  = f[i]
            l[i]  = l[i].replace("prime", "'")
            l[i]  = l[i].replace("p", ".")
        return n,p1,p2,l

    def mag_h(n):
        # Returns header and labels for color and magnitude plots
        if (n==1):
            lh=['Mbol','FUV_GALEX_AB','NUV_GALEX_AB','U_Johnson','B3_Johnson','V_Johnson','R_Cousins','I_Cousins','u_SDSS_AB','g_SDSS_AB','r_SDSS_AB',
                'i_SDSS_AB','z_SDSS_AB','K_Palomar','I5p7_IRAC','M24_MIPS']
            ly = ['M$_{BOL}$', 'GALEX FUV', 'GALEX NUV', 'U', 'B', 'V', 'R', 'I', 'u', 'g', 'r', 'i', 'z', 'K', 'IRAC 5.7 $\mu$m', 'MIPS 24 $\mu$m']
        elif (n==2):
            lh = ['J_2Mass', 'H_2Mass', 'Ks_2Mass', 'K_Palomar', 'L_Johnson', 'I3p6_IRAC', 'I4p5_IRAC', 'I5p7_IRAC', 'I7p9_IRAC', 'I12_IRAS', 'M24_MIPS',
                'I25_IRAS', 'I60_IRAS', 'M70_MIPS', 'I100_IRAS', 'M160_MIPS' ]
            ly = ['2Mass J', '2Mass H', '2Mass Ks', 'K', 'L', 'IRAC 3.6 $\mu$m', 'IRAC 4.5 $\mu$m', 'IRAC 5.7 $\mu$m', 'IRAC 7.9 $\mu$m', 'IRAS 12 $\mu$m',
                'MIPS 24 $\mu$m', 'IRAS 25 $\mu$m', 'IRAS 60 $\mu$m', 'MIPS 70 $\mu$m', 'IRAS 100 $\mu$m', 'MIPS 160 $\mu$m']
        return lh,ly

    print()
    s=n
    if (n==0):
        n,f1,l = umag_1(m1)
        if n > 0:
            nplt(n,a1,f1,l,file1,True,save)

    elif (n<=2):
        lh,ly = mag_h(n)
        splt(4,4,a1,m1,lh,ly,file1,True,save)

    elif (n==3):
        n,f1,f2,l = umag_2(m1,m2)
        if n > 0:
            nplt2(n,a1,f1,a2,f2,l,file1,file2,True,save)

    elif (n<=5):
        lh,ly = mag_h(n-3)
        splt2(4,4,a1,m1,a2,m2,lh,ly,file1,file2,True,save)

    if n > 0 and bt.sc:
        # Rename saved file
        if s <= 3:
            filout = os.environ.get('glxout') + '/' + bs.lnam(file1).replace('fits','mags') + str(s) + '.png'
        else:
            filout = os.environ.get('glxout') + '/' + bs.lnam(file1).replace('fits','') + bs.lnam(file2).replace('fits','mags') + str(s) + '.png'
        filtmp = os.environ.get('glxtmp') + '/' + 'tmpplot.png'
        os.rename(filtmp,filout)
        print(bs.gp(' Plot stored in file: ' + filout))
    return n

def normwav(o):
    # Returns normalization wavelength
    om[9]=''
    wn = 0.
    if o=='n' or o=='N':
        print()
        wx = input('Normalize SED''s at wavelength (A), press ENTER for no normalization = ')
#   elif o=='N':
#       print()
#       wx = input('Normalize SED''s at wavelength (A), press ENTER for no normalization = ')
    else:
        wx=''
    if len(wx) > 0 and float(wx) > 0:
        om[9] = wx
        wn = float(wx)
    return wn

def nplt(k,t,f,l,file1,flip,save):
    # Plot figure with n subplots

    def nvw(ax,i,k,t,f,l,flip):
        # Plot in window ax[i]
        if k > 1:
            bx = ax[i]
        else:
            bx = ax
        bx.set_xscale('log')
        bx.set_xlabel("t (yr)")
        bx.set_ylabel(l)
        if flip:
            bx.invert_yaxis()
        bx.plot(t,f,'r',linestyle='solid')

    title = 'Model 1: ' + bs.lnam(file1).replace('.fits','')
    if k==1:
        fig, ax = plt.subplots(1,k,figsize=(12,8))
        ax.set_title(title,loc='left',fontsize=10)
    else:
        fig, ax = plt.subplots(1,k,figsize=(5*k,5))
        ax[0].set_title(title,loc='left',fontsize=10)
    for i in range (k):
        nvw(ax,i,k,t,f[i],l[i],flip)
    fig.tight_layout()
    if bt.sc:
        filetmp = os.environ.get('glxtmp') + '/' + 'tmpplot.png'
        plt.savefig(filetmp,facecolor='w')
    plt.show()

def nplt2(k,t1,f1,t2,f2,l,file1,file2,flip,save):
    # Plot figure with n subplots, 2 models

    def nvw2(ax,i,k,t1,f1,t2,f2,l,flip):
        # Plot in window ax[i]
        if k > 1:
            bx = ax[i]
        else:
            bx = ax
        bx.set_xscale('log')
        bx.set_xlabel("t (yr)")
        bx.set_ylabel(l)
        if flip:
            bx.invert_yaxis()
        bx.plot(t1,f1,'b',linestyle='solid',label='1')
        bx.plot(t2,f2,'r',linestyle='solid',label='2')
        bx.legend(loc='best')

    titl1 = 'Model 1: ' + bs.lnam(file1).replace('.fits','')
    titl2 = 'Model 2: ' + bs.lnam(file2).replace('.fits','')
    if k==1:
        fig, ax = plt.subplots(1,k,figsize=(12,8))
        ax.set_title(titl1,loc='left',fontsize=10)
        ax.set_title(titl2,loc='right',fontsize=10)
    else:
        fig, ax = plt.subplots(1,k,figsize=(5*k,5))
        ax[0].set_title(titl1,loc='left',fontsize=10)
        ax[k-1].set_title(titl2,loc='left',fontsize=10)
    for i in range (k):
        nvw2(ax,i,k,t1,f1[i],t2,f2[i],l[i],flip)
    fig.tight_layout()
    if bt.sc:
        filetmp = os.environ.get('glxtmp') + '/' + 'tmpplot.png'
        plt.savefig(filetmp,facecolor='w')
    plt.show()

def phy(n):
    # Plot physical properties for option n

    def uphy_1(m):
        # Plot user selected properties
        global bx, one
        # Units
        u = ['(yr)','(photons/Mo)','(photons/Mo)','(photons/Mo)','','','','','', 'SN/yr', 'SN/yr', 'SN/yr', 'SN/yr', 'SN/yr',
            'PN/yr', 'number/Mo', 'number/Mo', 'number/Mo', 'Lo', '(Mo)','(Mo)', '(Mo)','(Mo)','(Mo/yr)','(Mo)','(Mo/Lo)',
            '(Mo/Lo)','(Mo/Lo)','(Mo/Lo)','(Mo/Lo)','(Mo/Lo)','(Nstars/yr)','(Nstars/yr/Lo)','(Mo)','','(M/Mo/yr)','(M/Mo/yr)']
        if not bt.jupy:
            # Ask for properties to plot, 1 model
            if one==0:
                print()
                wg = input('Select up to 4 properties to plot (<Enter> exits) N1 N2 N3 N4 = ')
                wg = wg.replace(",", " ")
                wg = wg.split()
                n = min(4,len(wg))
            else:
                wg = [str(one)]
                n = 1
        else:
            wg = om[7]
            n = min(4,len(wg))
        f = [""]*n
        for i in range(n):
            f[i] = bx[int(wg[i])]
        #print(i,a[i],f[i])
        if n <= 0:
            p=0.
            l=''
            n=-1
            return n,p,l
        # k bands to plot
        p = [0.]*n
        l = [""]*n
        for i in range (n):
            l[i] = f[i] + ' ' + u[int(wg[i])]
            l[i] = l[i].replace("_", "/")
            p[i] = m[f[i]]
        return n,p,l

    def uphy_2(m1,m2):
        # Plot user selected properties
        global bx, one
        # Units
        u = ['(yr)','(photons/Mo)','(photons/Mo)','(photons/Mo)','','','','','', 'SN/yr', 'SN/yr', 'SN/yr', 'SN/yr', 'SN/yr',
            'PN/yr', 'number/Mo', 'number/Mo', 'number/Mo', 'Lo', '(Mo)','(Mo)', '(Mo)','(Mo)','(Mo/yr)','(Mo)','(Mo/Lo)',
            '(Mo/Lo)','(Mo/Lo)','(Mo/Lo)','(Mo/Lo)','(Mo/Lo)','(Nstars/yr)','(Nstars/yr/Lo)','(Mo)','','(M/Mo/yr)','(M/Mo/yr)']
        if not bt.jupy:
            # Ask for properties to plot, 2 models
            if one==0:
                print()
                wg = input('Select up to 4 properties to plot (<Enter> exits) N1 N2 N3 N4 = ')
                wg = wg.replace(",", " ")
                wg = wg.split()
                n = min(4,len(wg))
            else:
                wg = [str(one)]
                n = 1
        else:
            wg = om[7]
            n = min(4,len(wg))
        f = [""]*n
        for i in range(n):
            f[i] = bx[int(wg[i])]
        if n <= 0:
            p1=0.
            p2=0.
            l=''
            n=-1
            return n,p1,p2,l
        # k bands to plot
        p1 = [0.]*n
        p2 = [0.]*n
        l  = [""]*n
        for i in range (n):
            l[i]  = f[i] + ' ' + u[int(wg[i])]
            l[i]  = l[i].replace("_", "/")
            p1[i] = m1[f[i]]
            p2[i] = m2[f[i]]
        return n,p1,p2,l

    def phy_h(n):
        # Returns header and labels for physical property plots
        if (n == 1):
            # Galaxy mass and M/L ratios
            lh = ['Mgalaxy', 'Mstars', 'Mremnants', 'Mretgas', 'Mtot', 'Mtot_Lb', 'Mtot_Lv', 'Mtot_Lk', 'SFR', 'Mliv_Lb', 'Mliv_Lv', 'Mliv_Lk', 'EvolutionaryFlux',
                'SpecificEvolFlux', 'TurnoffMass', 'BPMS_BMS']
            ly = ['Mgal (M$_\odot$)', 'Mliv (M$_\odot$)', 'Mrem (M$_\odot$)', 'Mgas (M$_\odot$)', 'Mtot = Mliv + Mrem (M$_\odot$)', 'Mtot/L$_B$ (M$_\odot$/L$_\odot$)',
                'Mtot/L$_V$ (M$_\odot$/L$_\odot$)', 'Mtot/L$_K$ (M$_\odot$/L$_\odot$)', 'SFR (M$_\odot$/yr)', 'Mliv/L$_B$ (M$_\odot$/L$_\odot$)', 'Mliv/L$_V$ (M$_\odot$/L$_\odot$)',
                'Mliv/L$_K$ (M$_\odot$/L$_\odot$)', 'Evolutionary Flux (Nstars/yr)', 'Specific Evol. Flux (Nstars/yr/L$_\odot$)', 'MS Turnoff Mass (M$_\odot$)',
                'Ratio PostMS/MS Bol. Flux']
        elif (n == 2):
            # Ionizing photons and other physical properties
            lh = ['logNLy', 'logNHeI', 'logNHeII', 'NBH', 'NNS', 'NWD', 'B912', 'SNR', 'PNBR']
            ly = ['log N(HI) (photons/M$_\odot$)', 'log N(HeI) (photons/M$_\odot$)', 'log N(HeII) (photons/M$_\odot$)', 'N(BH) (number/M$_\odot$)', 'N(NS) (number/M$_\odot$)',
                'N(WD) (number/M$_\odot$)', 'Lyman Break Amplitude', 'SN rate (SN/yr/L$_\odot$/M$_\odot$)', 'PNBR (PN/yr/L$_\odot$/M$_\odot$)' ]
        return lh,ly

    print()
    s=n
    if (n==0):
        n,f1,l = uphy_1(p1)
        if n > 0:
            nplt(n,a1,f1,l,file1,False,save)

    elif (n==1):
        lh,ly=phy_h(n)
        splt(4,4,a1,p1,lh,ly,file1,False,save)

    elif (n==2):
        lh,ly=phy_h(n)
        splt(3,3,a1,p1,lh,ly,file1,False,save)

    elif (n==3):
        n,f1,f2,l = uphy_2(p1,p2)
        if n > 0:
            nplt2(n,a1,f1,a2,f2,l,file1,file2,False,save)

    elif (n==4):
        lh,ly=phy_h(n-3)
        splt2(4,4,a1,p1,a2,p2,lh,ly,file1,file2,False,save)

    elif (n==5):
        lh,ly=phy_h(n-3)
        splt2(3,3,a1,p1,a2,p2,lh,ly,file1,file2,False,save)

    if n > 0 and bt.sc:
        # Rename saved file
        if s<=3:
            filout = os.environ.get('glxout') + '/' + bs.lnam(file1).replace('fits','physp') + str(s) + '.png'
        else:
            filout = os.environ.get('glxout') + '/' + bs.lnam(file1).replace('fits','') + bs.lnam(file2).replace('fits','physp') + str(s) + '.png'
        filtmp = os.environ.get('glxtmp') + '/' + 'tmpplot.png'
        os.rename(filtmp,filout)
        print(bs.gp(' Plot stored in file: ' + filout))
    return n

def sed_h(n):
    # Returns time header for age array in fits file
    t = 0.
    if n>=4:
        n=n-4
    if n==0:
        a = om[7]
        if ':' not in a:
            a = np.array(a, dtype=np.float32)
            if a[0] < 0:
                t=-a[0]
                if t<=0.1:
                    n=1
                elif t <=1:
                    n=2
                else:
                    n=3
        else:
            t  = 0.
            c  = a.replace(':',' ')
            c  = c.split()
            if '-' not in c[0]:
                c0 = float(c[0])
                c1 = float(c[1])
                c  = [c0,c1]
                k,i,b = geth(c,t1)
                k1 = int(i[1]) - 1
                k2 = int(i[2])
            else:
                k1 = -int(c[0])
                k2 =  int(c[1])
            a=[]
            for k in range(k1,k2):
                a.append(str(t1[k]*1.E-9))

    if n==1:
    	a = ['0.000015','0.00005','0.0001','0.0005','0.001','0.005','0.01','0.02','0.03','0.04','0.05','0.06','0.07','0.08','0.09','0.1']
    elif n==2:
    	a = ['0.001','0.005','0.01','0.025','0.05','0.075','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1']
    elif n==3:
    	a = ['0.01','0.05','0.1','0.25','0.5','0.75','1','2','3','4','5','6','7','8','9','10','11','12','13','14']
    a = np.array(a, dtype=np.float32)
    a = np.unique(a)  # unique sorts array a and suppresses duplicate entries (no need for a = np.sort(np.unique(a)))
    if t>0:
        a = a[a <= t]
    return a

def geth(a,t):
    # Find number of records h[k] corresponding to age a[k] (in Gyr) selected by user
    n = len(a)				# Number of time steps requested
    i = ['']*(n+1)
    b = ['']*(n+1)
    k = 0
    for l in range(n):
        k = k+1
        c = a[l]
        if c > 1.E3:
            c = c*1.E-9			# if age entered in yr transform to Gyr
        j = np.abs(t-c*1.E9).argmin()
        i[k] = int(j)+1
        if c >= 1:
            c = round(c*100)/100.
            c = f'{c:5.2f}'
            c = c.replace('.00','').replace('.50','.5').replace(' ','')
            c = c.replace('.90','.9').replace('.80','.8').replace('.70','.7').replace('.60','.6')
            c = c.replace('.40','.4').replace('.30','.3').replace('.20','.2').replace('.10','.1')
            b[k] = str(c)
        else:
            b[k] = str(a[l])
    #print(n,k,i,b)
    return k,i,b

def addext(filtmp,filout):
    # Adds unique extension to plot file name
    for i in range(100):
        if i < 10:
            ext = '.00' + str(i) + '.png'
        else:
            ext = '.0'  + str(i) + '.png'
        if not os.path.isfile(filout + ext):
            os.rename(filtmp,filout+ext)
            return filout+ext

def asplit(a):
    # Splits a into 2 lists according to character ';'
    if ';' not in str(a):
        if '|' in a:
            a1 = a.replace('|',' ')
            a1 = a1.split()
        else:
            if ':' in a:
                a1 = a
            else:
                a1 = [a]
        a2 = a1
    else:
        i  = a.find(';')
        a1 = a[:i]
        if '|' in a1:
            a1 = a1.replace('|',' ')
            a1 = a1.split()
        else:
            if ':' not in a1:
                a1 = [a1]
        a2 = a[i+1:]
        if '|' in a2:
            a2 = a2.replace('|',' ')
            a2 = a2.split()
        else:
            if ':' not in a2:
                a2 = [a2]
    return a1,a2

def sed(n):
    # Plots various sed's from the same model in a single panel

    # Detect if BC03 or CB20 model according to last time step
    s=n
    a1, a2 = asplit(om[7])
    om[7] = a1
    a  = sed_h(n)
    if n <= 3:
        # Get ages required
        e1 = round(t1[len(t1)-1]*1.E-9)	# Last time step rounded to nearest integer in Gyr
        e1 = max(0.2,e1)
        if e1 > 16:
            b1 = True                 # BC2003 models reach 20 Gyr
        else:
            b1 = False                # CB2020 models reach 14 Gyr
        r1 = a[a <= e1]
        # Get column header of sed's to plot
        k1,i1,b1 = geth(r1,t1)
    else:
        # Get ages required
        e1 = round(t1[len(t1)-1]*1.E-9)	# Last time step rounded to nearest integer in Gyr
        e2 = round(t2[len(t2)-1]*1.E-9)	# Last time step rounded to nearest integer in Gyr
        e1 = max(0.2,e1)
        e2 = max(0.2,e2)
        if e1 > 16:
            b1 = True                 # BC2003 models reach 20 Gyr
        else:
            b1 = False                # CB2020 models reach 14 Gyr
        if e2 > 16:
            b2 = True                 # BC2003 models reach 20 Gyr
        else:
            b2 = False                # CB2020 models reach 14 Gyr
        r1 = a[a <= e1]
        om[7] = a2
        a  = sed_h(n)
        r2 = a[a <= e2]
        # Get column header of sed's to plot
        k1,i1,b1 = geth(r1,t1)
        k2,i2,b2 = geth(r2,t2)

    # Limits entered by user
    if len(om[8]) > 0:
        xmin = float(om[8][0])
        xmax = float(om[8][1])
        if xmin < 0:
            # Fnu requested
            nu = True
            xmin = -xmin
        else:
            nu = False
        if len(om[8]) > 2:
            ymin = 10.**float(om[8][2])
            ymax = 10.**float(om[8][3])
            ido  = False
        else:
            ymin = 0.
            ymax = 0.
            ido  = True
    else:
        xmin, xmax, ymin, ymax = 0., 0., 0., 0.
        nu = False
        ido= False
    bt.nu = nu

    if s <= 3:
        # Check normalization
        if len(om[9]) > 0:
            wn = float(om[9])
            if wn > 0.:
                # Normalize sed's, find position in array w of the normalization wavelength = wn
                n = True
                ln = np.searchsorted(w1, wn)
        else:
            n = False
            wn = 0.

        # Single plot for model 1
        fig, ax = plt.subplots(1,1,figsize=(16,8))
        ax.minorticks_on()
        title = bs.lnam(file1).replace('.fits','') ; title.replace('.ised','')
        ax.set_title(title)
        ax.set_xlabel("Wavelength ($\AA$)")
        if n:
            if nu:
                ax.set_ylabel('L$\,_\\nu$ / L$\,_\\nu\,$(' + str(int(wn)) + '$\AA$)')
            else:
                ax.set_ylabel('L$\,_\lambda$ / L\,$_\lambda\,$(' + str(int(wn)) + '$\AA$)')
        else:
            if nu:
                ax.set_ylabel("L$\,_\\nu$ (L$_\odot$/Hz)")
            else:
                ax.set_ylabel("L$\,_\lambda$ (L$_\odot/\AA$)")
        plt.grid(True)
        if xmin != xmax:
            # Set plot limits entered by user
            if xmax > 1.E4:
                ax.set_xscale('log')
            ax.set_xlim(xmin,xmax)
            if ymin != ymax:
                ax.set_ylim(ymin,ymax)
            else:
                imin = np.searchsorted(w1, xmin)
                imax = np.searchsorted(w1, xmax)
                fmin = +1.E32
                fmax = -1.E32
        else:
            # Plot sed's in full range of wavelength and flux
            ax.set_xscale('log')
        ax.set_yscale('log')

        # Check if Fnu has been requested
        if nu:
            fc= 1.E-8/3.E10
            cn = fc*w1**2
        else:
            cn = [1.]*(len(w1))

        # Plot SED
        for l in range(1,k1+1):
            if n:
                # Normalize sed's
                if qf1:
                    fn = f1[ln][i1[l]]
                else:
                    fn = f1[i1[l]][ln]
                yp = f1[i1[l]][:]/fn/cn[ln]
            else:
                yp = f1[i1[l]][:]
            yp = cn*yp
            if ido:
                fmin= min(fmin,min(yp[imin:imax]))
                fmax= max(fmax,max(yp[imin:imax]))
            ax.plot(w1,yp,linestyle='solid',label=b1[l])
        if ido:
            ax.set_ylim(0.9*fmin,1.1*fmax)
        ax.legend(loc='best', title='Age (Gyr)', frameon=True)
        if bt.ss:
            # Save and rename file
            filtmp = os.environ.get('glxtmp') + '/' + 'tmpplot.png'
            plt.savefig(filtmp,facecolor='w')
            filout = os.environ.get('glxout') + '/' + bs.lnam(file1).replace('fits','sed')
            filout = addext(filtmp,filout)
        plt.show()
        if bt.ss:
            print(bs.gp(' Plot stored in file: ' + filout))

    else:
        # Plot sed's corresponding to each model in different frames
        # Check normalization
        if len(om[9]) > 0:
            wn = float(om[9])
            if wn > 0.:
                # Normalize sed's, find position in arrays w1,w2 of the normalization wavelength = wn
                n = True
                l1 = np.searchsorted(w1, wn)
                l2 = np.searchsorted(w2, wn)
        else:
            n = False
            wn = 0.

        # Check if Fnu has been requested
        if nu:
            fc= 1.E-8/3.E10
            c1 = fc*w1**2
            c2 = fc*w2**2
        else:
            c1 = [1.]*(len(w1))
            c2 = [1.]*(len(w2))

        # Find limits in common to both data groups
        fy1 = 1.E32 ; fy2 = -fy1
        if ido:
            imin = np.searchsorted(w1, xmin)
            imax = np.searchsorted(w1, xmax)
            jmin = np.searchsorted(w2, xmin)
            jmax = np.searchsorted(w2, xmax)
        for l in range(1,k1+1):
            if n:
                # Normalize sed's
                if qf1:
                    fn = f1[l1][i1[l]]
                else:
                    fn = f1[i1[l]][l1]
                yp = f1[i1[l]][:]/fn/c1[l1]
            else:
                yp = f1[i1[l]][:]
            yp = c1*yp
            if ido:
                fy1= min(fy1,min(yp[imin:imax]))
                fy2= max(fy2,max(yp[imin:imax]))
            else:
                fy1= min(fy1,min(yp))
                fy2= max(fy2,max(yp))
        for l in range(1,k2+1):
            if n:
                # Normalize sed's
                if qf2:
                    fn = f2[l2][i2[l]]
                else:
                    fn = f2[i2[l]][l2]
                yp = f2[i2[l]][:]/fn/c2[l2]
            else:
                yp = f2[i2[l]][:]
            yp = c2*yp
            if ido:
                fy1 = min(fy1,min(yp[jmin:jmax]))
                fy2 = max(fy2,max(yp[jmin:jmax]))
            else:
                fy1 = min(fy1,min(yp))
                fy2 = max(fy2,max(yp))
        fy1 = 0.9*max(fy1,1.e-25)
        fy2 = 1.1*fy2
        fx1 = 0.9*min(w1[0],w2[0])
        fx2 = 1.1*max(w1[len(w1)-1],w2[len(w2)-1])

        # Define plotting area
        if k1 > 1 or k2 > 1:
            fig, ax = plt.subplots(1,2,figsize=(20,8))
            ax0 = ax[0]
            ax1 = ax[1]
            one = False
            lbl = False
        else:
            fig, ax = plt.subplots(1,1,figsize=(16,8))
            ax0 = ax
            ax1 = ax
            one = True
            lbl = True

        # Plot first fits file on the LHS panel
        ax0.minorticks_on()
        ax0.set_xlabel('Wavelength ($\AA$)')
        if n:
            if nu:
                ax0.set_ylabel('L$\,_\\nu$ / L$\,_\\nu\,$(' + str(int(wn)) + '$\AA$)')
            else:
                ax0.set_ylabel('L$\,_\lambda$ / L$\,_\lambda\,$(' + str(int(wn)) + '$\AA$)')
        else:
            if nu:
                ax0.set_ylabel("L$\,_\\nu$ (L$_\odot$/Hz)")
            else:
                ax0.set_ylabel("L$\,_\lambda$ (L$_\odot/\AA$)")
        title = bs.lnam(file1).replace('.fits','') ; title.replace('.ised','')
        if one:
            title = '1: ' + title
            ax0.set_title(title,loc='left')
        else:
            ax0.set_title(title)
        ax0.set_yscale('log')
        if xmin != xmax:
            # Set plot limits entered by user
            if xmax > 1.E4:
                ax0.set_xscale('log')
            ax0.set_xlim(xmin,xmax)
            if ymin != ymax:
                ax0.set_ylim(ymin,ymax)
            else:
                ax0.set_ylim(fy1,fy2)
        else:
            # Plot sed's in full range of wavelength and flux
            if fx2 > 1.E4:
                ax0.set_xscale('log')
            ax0.set_xlim(fx1,fx2)
            ax0.set_ylim(fy1,fy2)
            ax0.set_xscale('log')
        ax0.set_yscale('log')
        for l in range(1,k1+1):
            if n:
                # Normalize sed's
                if qf1:
                    fn = f1[l1][i1[l]]
                else:
                    fn = f1[i1[l]][l1]
                yp = f1[i1[l]][:]/fn/c1[l1]
            else:
                yp = f1[i1[l]][:]
            yp = c1*yp
            if one:
               #ax0.plot(w1,yp,linestyle='solid',c='black',label=b1[l],zorder=1)
                ax0.plot(w1,yp,linestyle='solid',c='black',label=b1[l] + ' (model 1)',zorder=1)
            else:
                ax0.plot(w1,yp,linestyle='solid',label=b1[l])
        ax0.legend(loc='best', title='Age (Gyr)', frameon=True)
        xmin, xmax = ax0.set_xlim()  # return the current xlim
        ymin, ymax = ax0.set_ylim()  # return the current ylim
        ax0.grid(True)

        # Plot second fits file on the RHS panel
        ax1.minorticks_on()
        ax1.set_xlabel('Wavelength ($\AA$)')
        if n:
            if nu:
                ax1.set_ylabel('L$\,_\\nu$ / L$\,_\\nu\,$(' + str(int(wn)) + '$\AA$)')
            else:
                ax1.set_ylabel('L$\,_\lambda$ / L\,$_\lambda\,$(' + str(int(wn)) + '$\AA$)')
        else:
            if nu:
                ax1.set_ylabel("L$\,_\\nu$ (L$_\odot$/Hz)")
            else:
                ax1.set_ylabel("L$\,_\lambda$ (L$_\odot/\AA$)")
        title = bs.lnam(file2).replace('.fits','') ; title.replace('.ised','')
        if one:
            title = '2: ' + title
            ax1.set_title(title,loc='right')
        else:
            ax1.set_title(title)
        ax1.set_xlim(xmin,xmax)
        ax1.set_ylim(ymin,ymax)
        if xmax > 1.E4:
            ax1.set_xscale('log')
        ax1.set_yscale('log')
        for l in range(1,k2+1):
            if n:
                # Normalize sed's
                if qf2:
                    fn = f2[l2][i2[l]]
                else:
                    fn = f2[i2[l]][l2]
                yp = f2[i2[l]][:]/fn/c2[l2]
            else:
                yp = f2[i2[l]][:]
            yp = c2*yp
            if one:
               #ax1.plot(w2,yp,linestyle='dotted',c='red',label=b2[l],zorder=2)
                ax1.plot(w2,yp,linestyle='solid', c='red',label=b2[l] + ' (model 2)', zorder=2)
            else:
                ax1.plot(w2,yp,linestyle='solid',label=b2[l])
        ax1.legend(loc='best', title='Age (Gyr)', frameon=True)
        ax1.grid(True)
        if bt.ss:
            # Save and rename file
            filtmp = os.environ.get('glxtmp') + '/' + 'tmpplot.png'
            plt.savefig(filtmp,facecolor='w')
            filout = os.environ.get('glxout') + '/' + bs.lnam(file1).replace('fits','') + bs.lnam(file2).replace('fits','sed')
            filout = addext(filtmp,filout)
        plt.show()
        if bt.ss:
            print(bs.gp(' Plot stored in file: ' + filout))

def splt(i,j,t,m,lh,ly,file1,flip,save):
    # Plot figure with ixj subplots

    def vw(ax,i,j,t,m,k,lh,ly,flip):
        # Plots in window (i,j)
        if lh[k] =='':
            ax[i,j].set_axis_off()
        else:
            ax[i,j].set_xscale('log')
            ax[i,j].set_xlabel("t (yr)")
            ax[i,j].set_ylabel(ly[k])
            if flip:
                ax[i,j].invert_yaxis()
            if len(lh) > len(ly):
                ax[i,j].plot(t,m[lh[2*k]]-m[lh[2*k+1]],'r',linestyle='solid')
            else:
                ax[i,j].plot(t,m[lh[k]], 'r', linestyle='solid')

    fig, ax = plt.subplots(i,j,figsize=(10,8))
    k = -1
    for ii in range(i):
        for jj in range(j):
            ax[ii,jj].minorticks_on()
            k = k + 1
            vw(ax,ii,jj,t,m,k,lh,ly,flip)
    title = 'Model: ' + bs.lnam(file1).replace('.fits','')
    ax[0,0].set_title(title,loc='left',fontsize=10)
    for jj in range(j):
        fig.align_ylabels(ax[:,jj])       # align ylabel for the jj ax column
    fig.tight_layout()
    if bt.sc:
        filetmp = os.environ.get('glxtmp') + '/' + 'tmpplot.png'
        plt.savefig(filetmp,facecolor='w')
    plt.show()

def splt2(i,j,t1,m1,t2,m2,lh,ly,file1,file2,flip,save):
    # Plot figure with ixj subplots

    def vw2(ax,i,j,t1,m1,t2,m2,k,lh,ly,flip):
        # Plots in window (i,j)
        if lh[k] =='':
            ax[i,j].set_axis_off()
        else:
            ax[i,j].set_xscale('log')
            ax[i,j].set_xlabel("t (yr)")
            ax[i,j].set_ylabel(ly[k])
            if flip:
                ax[i,j].invert_yaxis()
            if len(lh) > len(ly):
                ax[i,j].plot(t1,m1[lh[2*k]]-m1[lh[2*k+1]],'b',linestyle='solid',label='1')
                ax[i,j].plot(t2,m2[lh[2*k]]-m2[lh[2*k+1]],'r',linestyle='solid',label='2')
            else:
                ax[i,j].plot(t1,m1[lh[k]],'b',linestyle='solid',label='1')
                ax[i,j].plot(t2,m2[lh[k]],'r',linestyle='solid',label='2')
        if i==0 and j==0:
            ax[i,j].legend(loc='best')
            #ax[i,j].legend(loc='best', title='Model:', frameon=True)

    fig, ax = plt.subplots(i,j,figsize=(10,8))
    k = -1
    for ii in range(i):
        for jj in range(j):
            ax[ii,jj].minorticks_on()
            k = k + 1
            vw2(ax,ii,jj,t1,m1,t2,m2,k,lh,ly,flip)
    for jj in range(j):
        fig.align_ylabels(ax[:,jj])       # align ylabel for the jj ax column
    # fig.suptitle("This Main Title is Nicely Formatted", fontsize=16)
    # fig.subplots_adjust(top=0.88)
    title = 'Model 1: ' + bs.lnam(file1).replace('.fits','')
    ax[0,0].set_title(title,loc='left',fontsize=10)
    title = 'Model 2: ' + bs.lnam(file2).replace('.fits','')
    ax[0,j-1].set_title(title,loc='right',fontsize=10)
   #title = 'Model 1: ' + bs.lnam(file1).replace('.fits','') + ',  Model 2: ' + bs.lnam(file2).replace('.fits','')
   #ax[0,0].set_title(title,loc='left',fontsize=8)
    fig.tight_layout()
    if bt.sc:
        filetmp = os.environ.get('glxtmp') + '/' + 'tmpplot.png'
        plt.savefig(filetmp,facecolor='w')
    plt.show()

def menu(k,n):
    # Plotting option menu
    if k==0:
        # CL version
        # Fits file options:
        v = ' '
        print()
        print('    Plot: s - SEDs                             Code: 0 - csp_galaxev             Legacy: 10 - model ID')
        print('          n - SEDs (normalized)                      1 - add_bursts                      11 - filter list')
        print('          m - photometric magnitudes                 2 - cm_evolution                    12 - ')
        print('          c - photometric colors                     3 - zmag                            13 - ')
        print('          i - line strength indices                  4 - galaxevpl')
        print('          p - physical properties                    5 - rf_phot')
        print('                                                     6 - of_phot')
        print('          '+bs.bp('h - help') + 22*v + '             7 - bcfits2txt')
        print('          l - list available models                  8 - ')
        print('          f - enter model file names                 9 - ')
        print('          q - quit')
        print()
        o = input(' choice = ')
    elif k<0:
        # CL version
        # ised and ascii file options:
        print()
        print(' Plot:    s - SEDs')
        print('          n - SEDs (normalized)')
        print('          f - enter new file name')
        print('          q - quit')
        print()
        o = input(' choice = ')
    else:
        # Jupyter Notebook version
        qp = om[2]
        op = om[3]
        if op == 'list models':
            o = 'l'
        elif op == 'Release Notes':
            o = 'h'
        elif op == 'filter list':
            o = '11'
        else:
            # Read required files
            n = rf(0,om)
            if op == 'csp_galaxev':
                o = '0'
            elif op == 'add_bursts':
                o = '1'
            elif op == 'cm_evolution':
                o = '2'
            elif op == 'zmag':
                o = '3'
            elif op == 'galaxevpl':
                o = '4'
            elif op == 'rf_phot':
                o = '5'
            elif op == 'of_phot':
                o = '6'
            elif op == 'bcfits2txt':
                o = '7'
            elif op == 'Model ID':
                o = '10'
            elif qp == 'SEDs':
                o = 's'
            elif qp == 'Normalized SEDs':
                o = 'n'
            elif qp== 'Magnitudes':
                o = 'm'
            elif qp == 'Colors':
                o = 'c'
            elif qp == 'Line indices':
                o = 'i'
            elif qp == 'Physical Properties':
                o = 'p'
            else:
                o = 's'
   #if o == 'q':
    if len(o) <= 0 or o == 'q':
       #quit()
       #sys.exit()
        return 0, 'q', 'q'
    elif o == 'l':
        bs.ldir('l')
    elif o == 'h':
        bs.myhelp()
        bs.ldir('?')
    elif o == '?':
        bs.myhelp()
        bs.ldir('?')
    elif o == 'f':
        # ask for file names
        n,om[0],om[1] = mfile()
        n = rf(0,om)
    elif o == '11':
        fl.filterid(0)
    q = qo(o)
    return n, o, q

def qo(o):
    # Assign parameter q for option o (indicates number of files allowed or required)
    if o == 's':
        q = 2
    elif o == 'n':
        q = 2
    elif o == 'm':
        q = 2
    elif o == 'c':
        q = 2
    elif o == 'i':
        q = 2
    elif o == 'p':
        q = 2
    elif o == 'l':
        q = 0
    elif o == 'f':
        q = 1
    elif o == 'h':
        q = 0
    elif o == '?':
        q = 0
    elif o == '0':
        q = 1
    elif o == '1':
        q = 1
    elif o == '2':
        q = 1
    elif o == '3':
        q = 1
    elif o == '4':
        q = 1
    elif o == '5':
        q = 1
    elif o == '6':
        q = 1
    elif o == '7':
        q = 1
    elif o == '10':
        q = 2
    elif o == '11':
        q = 0
    else:
        print('Unknown option: ',o)
        q = 'q'
    return q

def gg(o,m):
    # Displays first plot

    def hd(n):
        # Reports content of fits file
        print()
        print("In file " + bs.bp(bs.lnam(file1)) + " there are " + str('{:d}'.format(len(t1))) + " galaxy SED's, ranging from " + str('{:.1E}'.format(t1[0])) +
            " to " + str('{:.1E}'.format(t1[len(t1)-1])) + " yr")
        print("Each SED covers the wavelength range from " + str('{:.1f}'.format(w1[0])) + " to " + str('{:.1E}'.format(w1[len(w1)-1])) +
            " A in " + str('{:d}'.format(len(w1))) + " steps")
        if n==2:
            print()
            print("In file " + bs.bp(bs.lnam(file2)) + " there are " + str('{:d}'.format(len(t2))) + " galaxy SED's, ranging from " + str('{:.1E}'.format(t2[0])) +
                " to " + str('{:.1E}'.format(t2[len(t2)-1])) + " yr")
            print("Each SED covers the wavelength range from " + str('{:.1f}'.format(w2[0])) + " to " + str('{:.1E}'.format(w2[len(w2)-1])) +
            " A in " + str('{:d}'.format(len(w2))) + " steps")

    if m==0:
        return
    print()
    if o=='s':
        if m==1:
            sed(3)
        else:
            sed(7)
        file1 = bt.om[0]
        file2 = bt.om[1]
        hd(m)
    if o=='m':
        if m==1:
            mag(1)
        else:
            mag(4)
    if o=='c':
        if m==1:
            col(1)
        else:
            col(5)
    if o=='i':
        if m==1:
            idx(1)
        else:
            idx(7)
    if o=='p':
        if m==1:
            phy(1)
        else:
            phy(4)

def fplot(n,o,q):

    # Turn interactive mode on
    plt.ion()

    # Controls plotting routines
    if n==2 or q==1:
        # Plot selected options for one file
        n = fplot1(n,o)
    elif n==3 and q==2:
        # Plot selected options for two files
        n = fplot2(n,o)

def yl(l):
    # Builds y-label for variables in l
    n  = len(l)
    y = [" "]*(n+1)
    for i in range(n):
        lh = l[i]
        if lh=='':
            y[i] = lh
        elif ":m" in lh:
            lh = lh.replace(":m","")
            l[i] = lh
            y[i] = lh + " (mag)"
        elif ":n" in lh:
            lh = lh.replace(":n","")
            l[i] = lh
            y[i] = lh
        else:
            y[i] = lh + " ($\AA$)"
    return l,y

def mfile():
    # Ask for name of up to 2 fits file
    file1 = 'cb2019_z017_chab_hr_xmilesi_ssp.fits'
    n  = 1
    f1 = input(' BC_GALAXEV model 1 in file name [' + file1 + ']  = ')
    if f1 == 'q':
        quit()
    if len(f1) <= 0:
        f1 = bs.lnam(file1)
    else:
        f1 = bs.zrep(file1,f1)
    file1 = bs.fcheck(f1)
    file2 = bs.lnam(file1)
    f2    = input(' BC_GALAXEV model 2 in file name [' + file2 + ']  = ')
    if f1 == 'q':
        quit()
    if len(f2) <= 0:
        f2 = bs.lnam(f1)
    else:
        f2 = bs.zrep(f1,f2)
    file2 = bs.fcheck(f2)
    if file1 == file2:
        file2 = ''
        n=n+1
    else:
        n=3
    return n,file1,file2

def wcsp(file1):
    # Collects data to run csp_galaxev entered in widgetcsp
    # Input file name
    out = []
    file1 = bs.lnam(file1).strip()
    out.append(file1)
    # sfr
    io = rm[1] # ; o.write(io + '\n')
    to = rm[2] # ; o.write(to + '\n')
    # quenching
    if len(rm[5]) > 0:
        tcut = rm[5]
    else:
        tcut = '0'
    if io == '1':
        s = 'tau' + to + 'Gyr'
    elif io =='2':
        s = 'cmod' + to + 'Gyr'
    elif io =='3':
        s = 'cons' + to + 'Moyr'
    elif io =='4':
        io = '6'
        s = 'dlyd' + to + 'Gyr'
    elif io =='5':
        io = '8'
        s = 'linr' + to + 'Gyr'
    elif io =='6':
        io = '7'
        s = 'numr'
    out.append(io)
    out.append(to)
    out.append(tcut)
    out.append('n')
    # dust
    if len(rm[3]) > 0:
        out.append(True)
        out.append(rm[3][0])
        out.append(rm[3][1])
    else:
        out.append(False)
        out.append(0.)
        out.append(0.)
    # z rest frame
    if len(rm[4]) > 0:
        out.append(rm[4])
    else:
        out.append(0.)
    # output file
    if len(rm[6]) > 0:
        f = rm[6]
    else:
        f = bs.lnam(file1).replace('.fits','').replace('ssp',s).strip()
    out.append(f)
    print()
    print(bs.rp(' Running csp_galaxev code (may take a while)'))
    return out

def wadd():
    # Collects data to run add_bursts entered in widgetadd
    # Input file names,  burst parameters and output file
    file1 = bs.lnam(om[0]).strip()
    file2 = bs.lnam(om[1]).strip()
    print()
    print(bs.rp(' Running add_burst code (may take a while)'))
    b2  = rm[1].replace(",", " ").split()
    out = [ 0., 1., file1, float(b2[0])*1.E9, float(b2[1]), file2, rm[2] ]
    return out

def wcmev():
    # Collects data to run cm_evolution entered in widgetcmev
    # Input file name
    file1 = bs.lnam(om[0]).strip()
    # cosmology
    if len(rm[4]) > 0:
        a = rm[4] + ','
    else:
        a = '71' + ','
    if len(rm[5]) > 0:
        a = a + rm[5] + ','
    else:
        a = a + '0.27' + ','
    if len(rm[6]) > 0:
        a = a + rm[6] + ','
    else:
        a = a + '0.73'
    if len(rm[3]) > 0:
        b = rm[3]
    else:
        b = '13.5'
    # filters
    if len(rm[1]) > 0:
        c = rm[1]
    else:
        c = '15'
    if len(rm[2]) > 0:
        c = c + ',' + rm[2]
    else:
        c = c + ',' + '0'
    print()
    print(bs.rp(' Running cm_evolution code (may take a while)'))
    out = [ file1, a, b, c]
    return out

def wzmag(file1):
    # Collects data to run zmag entered in widgetzmag
    global f, t, fd, iread, h, q, clambda, tu, ttg, zf, omega, omega_lambda

    # Galaxy parameters read by widgets
    kf  = int(rm[1])
    kf  = [kf]
    z   = float(rm[2])
    ttg = float(rm[3])

    # Cosmological parameters entered in widgets
    h            = float(rm[4])
    omega        = float(rm[5])
    omega_lambda = float(rm[6])
    # Derive rest of parameters
    clambda,q    = cl.cosmol_c(h,omega,omega_lambda)
    tu           = cl.tuniverse(h,q,0.,clambda)
    zf           = cl.zx(ttg,h,q,clambda)
    dm           = cl.dismod(h,q,z)

    # Init variables
    bt.iread = True     # Read filter file only on first call
    p=-1                # Needed by 'percent' function on first call

    # Read BC/CB model in fits file
    w,f,t,e,*nouse = bc.bcfits(file1)
    hd.multhead(0,file1,t,w,e)

    # Compute redshift vs age
    evl.zage(t,ttg)

    # Get SED corresponding to redshift z = 0 (look-back-time = t0 = 0)
    t0, y0 = evl.zsed(0.,t,f)

    # Get SED corresponding to redshift z (look-back-time = tl)
    tl, yz = evl.zsed(z,t,f)
    tz=ttg-tl

    # Compute magnitudes vs. z
    evl.magz(w,z,kf,y0,yz,tl,tz,dm,file1)

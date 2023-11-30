import sys
import builtins      as bt
import basics.basics as b
import bcplot.bcplot as glx
import cspmod.cspmod as csp
import addmod.addmod as add
import bcevol.bcevol as evl
import bcfits.bcfits as bcf
import bcfilt.bcfilt as flt
import cosmol.cosmol as cml

def settings():
    # Read environment variables
    from dotenv import load_dotenv
    load_dotenv()
    # Define cosmological parameters
    bt.cosmol = cml.cosmol()
    #Init variables
    bt.lfits = ''       # Clear name of last fits file built by codes
    bt.fw    = 'cb2019_z017_chab_hr_xmilesi_ssp.fits'

def pyGALAXEV():
    settings()
    glx.pyGALAXEV()

def sspFileName():
    settings()
    glx.sspFileName()

def pyGALAXEVcl():
    settings()
    jupy = False
    save = False
    bt.save = save
    bt.jupy = jupy
    bt.ss   = False    # Don't store sed plots
    bt.sc   = False    # Don't store other plots

    if len(sys.argv) == 1:
        glx.pyGALAXEVcl(b.fr(['f','f'],0))

    elif sys.argv[1] == 'add' or sys.argv[1] == 'add_burts':
        add.pyadd_bursts(b.fr(sys.argv,0))

    elif sys.argv[1] == 'csp' or sys.argv[1] == 'csp_galaxev':
        csp.csp_galaxev(b.fr(sys.argv,0))

    elif sys.argv[1] == 'cmev' or sys.argv[1] == 'cm_evolution':
        evl.cmev(b.fr(sys.argv,0))

    elif sys.argv[1] == 'zmag':
        evl.zmag(b.fr(sys.argv,0))

    elif sys.argv[1] == 'rfphot' or sys.argv[1] == 'rf_phot':
        evl.krfp(b.fr(sys.argv,1))
        evl.rf_phot()

    elif sys.argv[1] == 'ofphot' or sys.argv[1] == 'of_phot':
        evl.kofp(b.fr(sys.argv,1))
        evl.of_phot()

    elif sys.argv[1] == 'gpl' or sys.argv[1] == 'galaxevpl':
        bcf.glxpl(b.fr(['f','f'],0))

    elif sys.argv[1] == 'addfilt' or sys.argv[1] == 'add_filters':
        flt.frm2fits()

    elif sys.argv[1] == 'bcf2t' or sys.argv[1] == 'bcfits2txt':
        bcf.bcf2t(b.fr(sys.argv,1))

    else:
        glx.pyGALAXEVcl(b.fr(sys.argv,1))

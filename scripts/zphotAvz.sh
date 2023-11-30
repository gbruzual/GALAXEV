#! /bin/zsh

./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av0.00  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av0.05  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av0.10  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av0.15  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av0.20  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av0.25  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av0.30  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av0.35  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av0.40  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av0.45  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av0.50  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av0.55  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av0.60  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av0.65  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av0.70  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av0.75  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av0.80  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av0.85  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av0.90  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av0.95  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av1.00  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av1.10  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av1.20  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av1.30  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av1.40  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av1.50  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av1.60  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av1.70  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av1.80  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av1.90  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av2.00  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av2.10  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av2.20  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av2.30  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av2.40  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av2.50  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av2.60  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av2.70  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av2.80  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av2.90  z$2  -o
./pyrfphot.py  cb2019_$1_chab_hr_xmilesi_ssp.fits  jpas  av3.00  z$2  -o
cd ../out
\rm -f cb2019_$1_chab_hr_xmilesi_ssp.jpas.z$2.ABmag cb2019_$1_chab_hr_xmilesi_ssp.jpas.z$2.Fjansky
cat cb2019_$1_chab_hr_xmilesi_ssp.jpas.z$2.Av*.ABmag   > cb2019_$1_chab_hr_xmilesi_ssp.jpas.z$2.ABmag
cat cb2019_$1_chab_hr_xmilesi_ssp.jpas.z$2.Av*.Fjansky > cb2019_$1_chab_hr_xmilesi_ssp.jpas.z$2.Fjansky
\rm cb2019_$1_chab_hr_xmilesi_ssp.jpas.z$2.Av*
cd ../pylib

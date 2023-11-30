#!/usr/bin/env python3

# Build fits table using astropy module

import sys
from astropy.table import Table

# Read ascii files and create fits file
f = sys.argv[1]
t = Table.read(f + '.sed', format='ascii')
t.write(f + '.fits', overwrite=True)  

t = Table.read(f + '.physical_properties',format='ascii')
t.write(f + '.fits', append=True)  

t = Table.read(f + '.photometry',format='ascii')
t.write(f + '.fits', append=True)  

t = Table.read(f + '.indices',format='ascii')
t.write(f + '.fits', append=True)  

t = Table.read(f + '.time_tsteps',format='ascii')
t.write(f + '.fits', append=True)  

t = Table.read(f + '.misc',format='ascii')
t.write(f + '.fits', append=True)  

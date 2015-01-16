#!/usr/bin/env python
#filename=regrid
#---------------------------------------by Gao John

#This program is designed to produce a new file that 
#contains regrided data from the original fits file.
#The regrid process is just simply take the average 
#value of every 4 pixels on the x-y plane in the original file.

from __future__ import division
from pylab import *
import pyfits
from sys import argv
import os
import time

#file name
fname = argv[1]
ind =fname.index('.',-5)
fname_o = fname[:ind]+'_rg.fits'

#assure there is no such a file
os.system('rm -f '+fname_o)

#open the file
hdulist = pyfits.open(fname)
co_header = hdulist[0].header
co_data = hdulist[0].data
[nn,ny,nx] = co_data.shape


#first cut image x y to make sure either can be divided by 2
print nx,ny
if (nx%2!=0): nx-=1
if (ny%2!=0): ny-=1
print nx,ny
#co_data = co_data[:,:ny-1,:nx-1]

new_data = zeros((nn,ny//2,nx//2))
for kk in range(nn):
    for jjj in range(ny//2):
        for iii in range(nx//2):
            jj = jjj*2
            ii = iii*2
            new_grid = [co_data[kk,jj,ii], co_data[kk,jj,ii+1], co_data[kk,jj+1,ii], co_data[kk,jj+1,ii+1]]
            new_data[kk,jjj,iii] = mean(new_grid)

hdulist[0].data = new_data
#synchronize header
if (hdulist[0].header['CRPIX1'] %2!=0):
    hdulist[0].header['CRVAL1'] -= 0.5*hdulist[0].header['CDELT1']
else:
    hdulist[0].header['CRVAL1'] += 0.5*hdulist[0].header['CDELT1']

hdulist[0].header['CRPIX1'] //= 2

if (hdulist[0].header['CRPIX2'] %2!=0):
    hdulist[0].header['CRVAL2'] -= 0.5*hdulist[0].header['CDELT2']
else:
    hdulist[0].header['CRVAL2'] += 0.5*hdulist[0].header['CDELT2']

hdulist[0].header['CRPIX2'] //= 2
hdulist[0].header['CDELT1'] *= 2.
hdulist[0].header['CDELT2'] *= 2.

#add history
time_now = time.asctime(time.localtime())
hdulist[0].header.add_history(time_now+'by python regrid of GY')
#write to
hdulist.writeto(fname_o)

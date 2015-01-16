#!/usr/bin/env python
#filename=cut_spec
#---------------------------------------by Gao John

#This program is designed to produce a new file that contains
#NO nan value point from a cube like fits file
#actully it may be only effective for CO spectra observing data

from __future__ import division
from pylab import *
import pyfits
from sys import argv
import os
import time

#as the name
def is_any_nan(kk):
    for ii in arange(nx):
        for jj in arange(ny):
            if(isnan(co_data[kk,jj,ii])):
                return True

    return False

#file name
fname = argv[1]
ind =fname.index('.',-5)
fname_o = fname[:ind]+'eff.fits'

#assure there is no such a file
os.system('rm -f '+fname_o)

#open the file
hdulist = pyfits.open(fname)
co_header = hdulist[0].header
co_data = hdulist[0].data

[nn,ny,nx] = co_data.shape
vn = ((arange(co_header['NAXIS3'])-co_header['CRPIX3'])*co_header['CDELT3']+co_header['CRVAL3'])/1000.

ocenter_x = co_header['CRVAL1']
ocenter_y = co_header['CRVAL2']
ocenter_v = co_header['CRVAL3']
ocenter_i = co_header['CRPIX1']
ocenter_j = co_header['CRPIX2']
ocenter_k = co_header['CRPIX2']
cdelt_x = co_header['CDELT1']
cdelt_y = co_header['CDELT2']
cdelt_v = co_header['CDELT3']

#first cut image x y
x_value = (arange(co_header['NAXIS1'])-ocenter_i)*cdelt_x+ocenter_x
y_value = (arange(co_header['NAXIS2'])-ocenter_j)*cdelt_y+ocenter_y
#============================================
raw_s = raw_input('Input the number of pixels in every directions, like:\'left_cut,right_cut,up_cut,bottom_cut\'')
l_cut = 0
r_cut = 0
u_cut = 0
b_cut = 0
l_cut,r_cut,u_cut,b_cut = map(int,raw_s.split(','))
#============================================

hdulist[0].data = co_data[:,b_cut:nx-u_cut,l_cut:ny-r_cut]#1:3->1,2
co_data = hdulist[0].data
[nn,ny,nx] = co_data.shape

#synchronize header
hdulist[0].header['CRPIX1'] = (l_cut+nx-1-r_cut)//2-l_cut
hdulist[0].header['CRPIX2'] = (b_cut+ny-1-u_cut)//2-b_cut
hdulist[0].header['CRVAL1'] = x_value[((l_cut+nx-1-r_cut)//2)]
hdulist[0].header['CRVAL2'] = y_value[(b_cut+ny-1-u_cut)//2]

kcut1 = 0
kcut2 = co_header['NAXIS3']-1
k_1or2 = 0

for kk in arange(nn):
    if(k_1or2==0):
        if(is_any_nan(kk)):
            kcut1 = kk+1
            continue

    k_1or2 = 1
    if(is_any_nan(kk)):
        kcut2 = kk-1
        break

print kcut1,kcut2
if (kcut1>=kcut2):
    print 'please check the value you just typed!'
    exit()

hdulist[0].data = co_data[kcut1:kcut2,:,:]

#synchronize header
hdulist[0].header['CRPIX3'] = (kcut1+kcut2)//2-kcut1
hdulist[0].header['CRVAL3'] = vn[(kcut1+kcut2)//2]*1000.

#add history
time_now = time.asctime(time.localtime())
hdulist[0].header.add_history(time_now+'by python cut_spec of GY')
#write to
hdulist.writeto(fname_o)

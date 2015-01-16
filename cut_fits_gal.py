#!/usr/bin/env python
#filename=cut_fits_gal
#---------------------------------------by Gaoyuan John
#a new version

#import sometings
from __future__ import division
from pylab import *
from sys import argv
import pyfits
import os

#input the filename that need to be cut

fname = argv[1]
ind = fname.index('.',-5)
fname_s = fname[:ind]+'_smaller.fits'

#assure there is no such a file
os.system('rm -f '+fname_s)

#----------------------------------------#

hdulist = pyfits.open(fname)
header = hdulist[0].header
data = hdulist[0].data

#galactic
#unit degree
center_x = 337.8
center_y = -0.1

#range need as a square
#unit degree
cut_range = 0.7

#----------------------------------------#

#what we need, 2pixel
#index [ii, jj]
ocenter_x = header['CRVAL1']
ocenter_y = header['CRVAL2']
ocenter_i = header['CRPIX1']
ocenter_j = header['CRPIX2']
cdelt_x = header['CDELT1']
cdelt_y = header['CDELT2']

ii = (arange(header['NAXIS1'])-ocenter_i)*cdelt_x+ocenter_x
jj = (arange(header['NAXIS2'])-ocenter_j)*cdelt_y+ocenter_y

x_top = center_x-cdelt_x/abs(cdelt_x)*cut_range/2.0
y_top = center_y-cdelt_y/abs(cdelt_y)*cut_range/2.
x_bottom = center_x+cdelt_x/abs(cdelt_x)*cut_range/2.
y_bottom = center_y+cdelt_y/abs(cdelt_y)*cut_range/2.

print x_top, x_bottom, y_top, y_bottom
a24 = 0
for i in range(header['NAXIS1']-2):
    if ((x_top>=ii[i] and x_top<=ii[i+1]) or (x_top<=ii[i] and x_top>=ii[i+1])):
      i_top = i
      a24+=1
    if ((y_top>=jj[i] and y_top<=jj[i+1]) or (y_top<=jj[i] and y_top>=jj[i+1])):
      j_top = i
      a24+=1
    if ((x_bottom>=ii[i] and x_bottom<=ii[i+1]) or (x_bottom<=ii[i] and x_bottom>=ii[i+1])):
      i_bottom = i
      a24+=1
    if ((y_bottom>=jj[i] and y_bottom<=jj[i+1]) or (y_bottom<=jj[i] and y_bottom>=jj[i+1])):
      j_bottom = i
      a24+=1

if (a24 != 4):
    print 'The boundaries collection failed!', ' a24=', a24, i_top, j_top, i_bottom, j_bottom
    exit(1)

if (i_top>=i_bottom or j_top>=j_bottom):
    print 'i:[', i_top, i_bottom, '], j:[', j_top, j_bottom, ']'
    exit(2)
else:
    print 'i:[', i_top, i_bottom, '], j:[', j_top, j_bottom, ']'

#----------------------------------------#
#nhdulist = pyfits.open('..//most//MOS348p5_smaller.FIT', mode = 'ostream')
#the order is reversed!!!!!!!!!!!!!
#print data[242:405][160:324]
#change data
hdulist[0].data = data[j_top:j_bottom+1, i_top:i_bottom+1]

#print hdulist[0].data
#change header
#nhdulist[0].header = header
hdulist[0].header['CRPIX1'] = (i_top+i_bottom)//2-i_top
hdulist[0].header['CRPIX2'] = (j_top+j_bottom)//2-j_top
hdulist[0].header['DATAMAX'] = hdulist[0].data.max()
hdulist[0].header['DATAMIN'] = hdulist[0].data.min()
#print dir(hdulist[0].data)

#hdulist[0].header['NAXIS1'] = i_bottom-i_top+1
#hdulist[0].header['NAXIS2'] = j_bottom-j_top+1

hdulist[0].header['CRVAL1'] = ii[(i_top+i_bottom)//2]
hdulist[0].header['CRVAL2'] = jj[(j_top+j_bottom)//2]

#print header
#----------------------------------------#

#print hdulist[0].header
#save file as
hdulist.writeto(fname_s)

#----------------------------------------#
#hdulist.flush()
#hdulist.close()

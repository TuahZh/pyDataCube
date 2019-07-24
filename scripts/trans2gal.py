#!/usr/bin/env python
#filename=trans2gal
#---------------------------------------by Gao John

from kapteyn import wcs
from numpy import *
import pyfits
from kapteyn.interpolation import map_coordinates
from sys import argv

#filename
fname = argv[1]
ind = fname.index('.',-5)
fname_o = fname[:ind]+'_gal.fits'

#open file ctb37b_12coall_MEAN.fits (12CO)
filein = pyfits.open(fname)
header = filein[0].header
#print header                                                   #check

proj1 = wcs.Projection(header)                                 #source projection
print proj1.skysys
trans = wcs.Transformation(proj1.skysys, skyout = wcs.galactic)

header['CTYPE1'], header['CTYPE2'] = 'GLON-TAN', 'GLAT-TAN'    #new axis types
header['CRVAL1'], header['CRVAL2'] = trans((header['CRVAL1'], header['CRVAL2']))                                                           
                                                               #new reference point

proj2 = wcs.Projection(header)                                 #destination

coords = wcs.coordmap(proj1, proj2)

image_in = filein[0].data
image_out = map_coordinates(image_in, coords, order=1, cval=NaN)

filein[0].data = image_out
filein.writeto(fname_o)

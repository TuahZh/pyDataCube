#!/usr/bin/env python
#filename=trans2fk5
#---------------------------------------by Gao John
#edit a new version

from kapteyn import wcs
from pylab import *
from sys import argv
import pyfits
from kapteyn.interpolation import map_coordinates
#
#filename
fname = argv[1]
ind = fname.index('.',-5)
fname_o = fname[:ind]+'_fk5.fits'

#
filein = pyfits.open(argv[1])
header = filein[0].header
#print header                                                   #check

proj1 = wcs.Projection(header)                                 #source projection
print proj1.skysys
trans = wcs.Transformation(proj1.skysys, skyout = (0, 6, 'J2000.0'))  #fk5, j2000

header['CTYPE1'], header['CTYPE2'] = 'RA---SIN', 'DEC--SIN'    #new axis types
header['CRVAL1'], header['CRVAL2'] = trans((header['CRVAL1'], header['CRVAL2']))                                                           #new reference point 

#correct header: epoch
del header['EPOCH']
header.update('EQUINOX', 2000, 'Equinox of equatorial coordinates', before='ORIGIN')
header.update('RADESYS', 'FK5', 'Equatorial coordinate system', before='EQUINOX')

proj2 = wcs.Projection(header)                                 #destination 
coords = wcs.coordmap(proj1, proj2) 

image_in = filein[0].data
image_out = map_coordinates(image_in, coords, order=3, cval=NaN) 
filein[0].data = image_out
filein.writeto(fname_o)

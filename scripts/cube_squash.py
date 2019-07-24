#########################################################################
# File Name: cube.py
# Author: ZHANG, GaoYuan
# mail: zgy0106@gmail.com
# Created Time: Mon 20 Jan 2014 10:13:17 PM CST
#########################################################################
#discription

#!/bin/env python

from __future__ import division
from pylab import *
import pyfits
import pywcsgrid2
from kapteyn import wcs
from mpl_toolkits.axes_grid1.axes_rgb import imshow_rgb
import gaoyuan

file_12co = pyfits.open("../kes41_12coall_MEANeff_based.fits")
file_13co = pyfits.open("../kes41_13coall_MEANeff_based.fits")
file_c18o = pyfits.open("../kes41_c18oall_MEANeff_based.fits")
cube12 = gaoyuan.Cube(file_12co[0].data, file_12co[0].header)
cube13 = gaoyuan.Cube(file_13co[0].data, file_13co[0].header)
cube18 = gaoyuan.Cube(file_c18o[0].data, file_c18o[0].header)

#=========================
v_min = -70e3
v_max = -40e3
#=========================

#squash
(da12, he12) = cube12.squash(v_min, v_max, coord="sky", mode="sum")
(da13, he13) = cube13.squash(v_min, v_max, coord="sky", mode="sum")
(da18, he18) = cube18.squash(v_min, v_max, coord="sky", mode="sum")

file_out12 = pyfits.PrimaryHDU(data=da12, header=he12)
file_out13 = pyfits.PrimaryHDU(data=da13, header=he13)
file_out18 = pyfits.PrimaryHDU(data=da18, header=he18)

file_out12.writeto("kes41_12_%d-%d.fits" % (v_min/1e3, v_max/1e3))
file_out13.writeto("kes41_13_%d-%d.fits" % (v_min/1e3, v_max/1e3))
file_out18.writeto("kes41_18_%d-%d.fits" % (v_min/1e3, v_max/1e3))

#draw
#
#ax = pywcsgrid2.subplot(111, header=file_12co[0].header)
#im = imshow_rgb(ax,
#        mytrim(r,0,100),
#        mytrim(g,0,100),
#        mytrim(b,0,100),
#        origin="lower")
#
#
#show()

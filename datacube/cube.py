#########################################################################
# File Name: gaoyuan.py
# Author: ZHANG, GaoYuan
# mail: zgy0106@gmail.com
# Created Time: Mon 20 Jan 2014 10:13:17 PM CST
#########################################################################
#discription

#!/bin/env python

#from __future__ import division
import numpy as np
try:
    from astropy.io import fits as pyfits
except:
    import pyfits
#import pywcsgrid2
from astropy.wcs import WCS as wcs
#
#========Class Cube===================
class Cube:
    def __init__(self, data, header):
        self.data = data
        self.header = header
        if (header["NAXIS"]==3):
            self.cube = data
        elif (header["NAXIS"]==4):
            self.cube = np.zeros(data.shape[1:3])
            self.cube = data[0,:,:,:]
        else:
            raise FITSTypeError


    def squash(self, index_min, index_max, coord="sky", mode="sum"):
        if (coord=="sky"):
            proj = wcs(self.header)
            line_ax = 3
            line = proj.sub(line_ax)
            slice_min = int(line.wcs_world2pix(index_min))
            slice_max = int(line.wcs_world2pix(index_max))

        if (mode=="sum"):
            (nn, ny, nx) = self.cube.shape
            sum = np.zeros((ny, nx))
            for i in xrange(nx):
                for j in xrange(ny):
                    for k in range(nn)[slice_min:slice_max]:
                        if (np.isnan(self.cube[k,j,i])):
                            continue

                        sum[j,i] += \
                            self.cube[k,j,i]*self.header["CDELT3"]/1000.

            header_new = self.header
            header_new.__delitem__("NAXIS4")
            header_new.__delitem__("CTYPE4")
            header_new.__delitem__("CRPIX4")
            header_new.__delitem__("CDELT4")
            header_new.__delitem__("CRVAL4")
            header_new["NAXIS"] = 2
            header_new.__delitem__("NAXIS3")
            header_new.__delitem__("CTYPE3")
            header_new.__delitem__("CRPIX3")
            header_new.__delitem__("CDELT3")
            header_new.__delitem__("CRVAL3")
            header_new.__delitem__("RESTFREQ")
            header_new.__delitem__("SPECSYS")
            header_new.__delitem__("DATAMAX")
            header_new.__delitem__("DATAMIN")
            header_new.add_history("squashed by GY")

        return sum, header_new

#    def channel_maps():


#=================Class Cube======================

#=================Function mytrim=================
def mytrim(d, vmin, vmax):
    dd = (d+vmin)/(vmin+vmax)
    return np.clip(dd, 0, 1)
#=================Function mytrim=================



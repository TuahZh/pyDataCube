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
        """Initialize a Cube
        cube is data
        if 3D data
        cube take the last 3 dimensions if 4D
        2D is not acceptable
        """
        self.data = data
        self.header = header
        if (header["NAXIS"]==3):
            self.cube = data
        elif (header["NAXIS"]==4):
            self.cube = np.zeros(data.shape[1:3])
            self.cube = data[0,:,:,:]
        else:
            raise FITSTypeError


    def squash(self, index_min, index_max, coord="sky", mode="sum", unit="km/s"):
        """
        This function is used to squash a datacube between a velocity/spectral
        interval.

        Currently only velocity is used and the original fits file is assumed to
        have a velocity unit of m/s

        Output:
            data     2D image
            header   Try to drop all the extra header item
        """
        index_min = self._unit_conv(index_min, unit=unit)
        index_max = self._unit_conv(index_max, unit=unit)

        if (coord=="sky"):
            (slice_min, slice_max) = self._vel2pix(index_min, index_max)

        if (mode=="sum"):
            (nn, ny, nx) = self.cube.shape
            sum = np.zeros((ny, nx))
            # This can be optimized
            for i in range(nx):
                for j in range(ny):
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
            header_new.__delitem__("CROTA4")
            header_new["NAXIS"] = 2
            header_new.__delitem__("NAXIS3")
            header_new.__delitem__("CTYPE3")
            header_new.__delitem__("CRPIX3")
            header_new.__delitem__("CDELT3")
            header_new.__delitem__("CRVAL3")
            header_new.__delitem__("CROTA3")
            header_new.__delitem__("RESTFREQ")
            header_new.__delitem__("SPECSYS")
            header_new.__delitem__("DATAMAX")
            header_new.__delitem__("DATAMIN")
            header_new.add_history("Squashed by pyDataCube")

        return sum, header_new

    def trimnan(self):
        pass

    def trim_val(self, index_min, index_max, coord="sky", unit="km/s"):
        index_min = self._unit_conv(index_min, unit=unit)
        index_max = self._unit_conv(index_max, unit=unit)

        if (coord=="sky"):
            (slice_min, slice_max) = self._vel2pix(index_min, index_max)
        pass

    def _unit_conv(val, unit="m/s"):
        if (unit=="km/s"):
            return val*1000. #unit as m/s
        elif (unit=="cm/s"):
            return val/100.
        elif (unit=="m/s"):
            pass
        else:
            raise ValueError("Cube cannot read this unit")

    def _vel2pix(self, vel_min, vel_max):
        proj = wcs(self.header)
        line_ax = ["spectral"]
        # sub is different from the original one
        line = proj.sub(line_ax)
        slice_min = line.wcs_world2pix(vel_min, 1)
        slice_max = line.wcs_world2pix(vel_max, 1)
        try:
            slice_min = int(slice_min[0].round())
            slice_max = int(slice_max[0].round())
        except:
            print("The version of astropy is not satisfied!")

        return (slice_min, slice_max)


#    def channel_maps():


#=================Class Cube======================

#=================Function mytrim=================
def mytrim(d, vmin, vmax):
    dd = (d+vmin)/(vmin+vmax)
    return np.clip(dd, 0, 1)
#=================Function mytrim=================



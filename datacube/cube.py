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
import warnings
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
        # link 
        self.data = data
        self.header = header
        # hard copy
        self._header_copy = header.copy()
        if (header["NAXIS"]==3):
            self.cube = data.copy()
        elif (header["NAXIS"]==4):
#            self.cube = np.zeros(data.shape[1:3])
            self.cube = (data[0,:,:,:]).copy()
        else:
            raise FITSTypeError


    def squash(self, index_min, index_max, coord="sky", mode="sum", unit="km/s"):
        """
        This function is used to squash a datacube between a velocity/spectral
        interval.

        Currently only velocity is used and the original fits file is assumed to
        have a velocity unit of m/s

        Include index_min and index_max (or velocity min and max)

        Output:
            data     2D image; unit is K*km/s
            header   Try to drop all the extra header item
        """
        index_min = self._unit_conv(index_min, unit=unit)
        index_max = self._unit_conv(index_max, unit=unit)

        if (coord=="sky"):
            (slice_min, slice_max) = self._vel2pix(index_min, index_max)

        if (mode=="sum"):
            (nn, ny, nx) = self.cube.shape
            sum_cube = np.zeros((ny, nx))
            # This can be optimized
#            for i in range(nx):
#                for j in range(ny):
#                    for k in range(nn)[slice_min:slice_max]:
#                        if (np.isnan(self.cube[k,j,i])):
#                            continue
#
#                        sum[j,i] += \
#                            self.cube[k,j,i]*self.header["CDELT3"]/1000.
            # include the final pix
            sum_cube = self.cube[slice_min:slice_max+1, : , :].sum(axis=0)* \
                self.header["CDELT3"]/1000. # The velocity has a unit m/s

            header_new = self.header.copy()
            try:
                header_new.__delitem__("NAXIS4")
                header_new.__delitem__("CTYPE4")
                header_new.__delitem__("CRPIX4")
                header_new.__delitem__("CDELT4")
                header_new.__delitem__("CRVAL4")
                header_new.__delitem__("CROTA4")
            finally:
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

        return sum_cube, header_new

    def trimnan(self):
        pass

    def trim_val(self, index_min, index_max, coord="sky", unit="km/s", inplace=False):
        """
        To trim the datacube according the interval of spectral velocity range
        or the channel number range

        Output:
            new_cube     new datacube with the spectral trim
            header       new header remain almost the same
        """
        index_min = self._unit_conv(index_min, unit=unit)
        index_max = self._unit_conv(index_max, unit=unit)

        if (coord=="sky"):
            (slice_min, slice_max) = self._vel2pix(index_min, index_max)

        #index_list = arange(self.cube.shape[0]) + 1. # fits has an origin from 1
        #cent_floor = np.floor(self.header["CRPIX3"])
        #cent_ceiling = np.ceiling(self.header["CRPIX3"])
        #index_floor = index_list[cent_floor]
        #index_ceiling = index_list[cent_ceiling]
        ### include the final pixel
        if (inplace):
            new_cube = self.cube[slice_min:slice_max+1, :, :]
            new_header = self.header
        else:
            new_cube = self.cube.copy()[slice_min:slice_max+1, :,:]
            new_header = self.header.copy()

        new_header.add_history("Trimed by pyDataCube")
        #new_header["CRVAL3"]
        #new_header["CDELT3"]
        #The central pixel should be changed
        with warnins.catch_warnins():
            warnings.filterwarnings("error")
            try:
                # slice_min and slice_max was for python, has a 0-origin
                new_header["CRPIX3"] = self._pix_trim(self.header["CRPIX3"], \
                                                      slice_min+1, slice_max+1, \
                                                      origin=1)
            except UserWarning:
                warnings.warn("Warning: changed central pix is not within the interval!")
                pass
                # here need something to do with CRVAL and CRPIX

        return new_cube, new_header

    def _unit_conv(self, val, unit="m/s"):
        if (unit=="km/s"):
            return val*1000. #unit as m/s
        elif (unit=="cm/s"):
            return val/100.
        elif (unit=="m/s"):
            pass
        else:
            raise RuntimeError("Cube cannot read this unit")

    def _vel2pix(self, vel_min, vel_max):
        proj = wcs(self.header)
        line_ax = ["spectral"]
        # sub is different from the original one
        line = proj.sub(line_ax)
        # slice is for numpy array should use 0-origin
        slice_min = line.wcs_world2pix(vel_min, 0)
        slice_max = line.wcs_world2pix(vel_max, 0)
        try:
            slice_min = int(slice_min[0].round())
            slice_max = int(slice_max[0].round())
        except:
            print("The version of astropy is not satisfied!")
            raise

        return (slice_min, slice_max)

#    def _pix_trim(self, d, vmin, vmax):
#        dd = (d+vmin)/(vmin+vmax)
#        return np.clip(dd, 0, 1)
     def _pix_trim(self, d, vmin, vmax, origin=1):
         """_pix_trim
         return a new index with the cut off edge vmin and vmax
         d must be in beween vmin and vmax? no
         d can be out of the range, but a warning is printed
         origin=1 means for Fits or Fortran
         origin=0 means for Python or C
         """
         if (np.isscaler(d)):
            if (d<vmin or d>vmax):
                warnings.warn("Warning: changed central pix is not within the interval!")
         else:
             if (np.logical_or(d<vmin,d>vmax).any()):
                warnings.warn("Warning: changed central pix is not within the interval!")

         dd = d-vmin+origin
         return dd

#    def channel_maps():


#=================Class Cube======================

#=================Function mytrim=================
#=================Function mytrim=================



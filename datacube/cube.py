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
from scipy.optimize import curve_fit
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
            except KeyError:
                pass
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

    def channel_squash(self, velo_min, velo_max, interval=1, unit="km/s", num=None):
        """
        Squash the data cube every certain interval
        currently only support km/s
        """
        if (num is not None):
            interval = (velo_max-velo_min)/num
        else:
            num = (velo_max-velo_min)/interval
            num = int(np.ceil(num))

        velo_mins = np.arange(num)*interval+velo_min
        velo_maxs = velo_mins+interval
        datas = np.zeros((num, self.cube.shape[1], self.cube.shape[2]))
        headers = {}
        for ii in range(num):
            datas[ii,:,:], headers[ii] = self.squash(velo_mins[ii], velo_maxs[ii])

        return datas, headers

    def trimnan(self, iterative=True, inplace=False, max_iter=5):
        """
        To trim the outer nan part from the datacube
        A ring of pixel is dropped every time until there is no nan left
        The 3rd dimension (spectral) is assumed to be all good without nan
        """
        if (inplace):
            new_cube = self.cube
            new_header = self.header
        else:
            new_cube = self.cube.copy()
            new_header = self.header.copy()

        num_iter = 0
        while (np.isnan(new_cube).any()):
            new_cube, new_header = self._trim_one(new_cube, new_header)
            num_iter += 1
            if (num_iter>max_iter):
                warnings.warn("The new file may still contain nan values!")
                break

        new_header.add_history("Trim the nan part in the datacube (by pyDataCube)")
        if (new_header["NAXIS"]==4):
            try:
                new_header.__delitem__("NAXIS4")
                new_header.__delitem__("CTYPE4")
                new_header.__delitem__("CRPIX4")
                new_header.__delitem__("CDELT4")
                new_header.__delitem__("CRVAL4")
                new_header.__delitem__("CROTA4")
            except KeyError:
                pass
            finally:
                new_header["NAXIS"] = 3

        return new_cube, new_header

    def _trim_one(self, ncube, nheader):
        """
        To drop the outmost ring of the cube
        """
        ny = ncube.shape[1]
        nx = ncube.shape[2]
        ncube = ncube[:, 1:-1, 1:-1]
        with warnings.catch_warnings():
            warnings.filterwarnings("error")
            try:
                # slice_min and slice_max was for python, has a 0-origin
                # the slice_max may not be right but not important
                nheader["CRPIX1"] = self._pix_trim(nheader["CRPIX1"], \
                                                      1+1, nx-1, \
                                                      origin=1)
                nheader["CRPIX2"] = self._pix_trim(nheader["CRPIX2"], \
                                                      1+1, ny-1, \
                                                      origin=1)
            except UserWarning:
                warnings.warn("Warning: changed central pix is not within the interval!")
                pass
                # here need something to do with CRVAL and CRPIX

        return ncube, nheader


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
        with warnings.catch_warnings():
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

        if (new_header["NAXIS"]==4):
            try:
                new_header.__delitem__("NAXIS4")
                new_header.__delitem__("CTYPE4")
                new_header.__delitem__("CRPIX4")
                new_header.__delitem__("CDELT4")
                new_header.__delitem__("CRVAL4")
                new_header.__delitem__("CROTA4")
            except KeyError:
                pass
            finally:
                new_header["NAXIS"] = 3


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
        if (np.isscalar(d)):
           if (d<vmin or d>vmax):
               warnings.warn("Warning: changed central pix is not within the interval!")
        else:
            if (np.logical_or(d<vmin,d>vmax).any()):
               warnings.warn("Warning: changed central pix is not within the interval!")

        dd = d-vmin+origin
        return dd

    def _get_line(self, pos, ncube=None, origin=1):
        """
        Get the spectral line from a position
        """
        if (ncube is None):
            ncube = self.cube

        if (len(pos)!=2):
            raise ValueError("Length of position should be 2!")

        self.catch_line = ncube[:, pos[1]-origin, pos[0]-origin]
        return self.catch_line

    def _get_average_line(self, ncube=None):
        """
        Get the average spectrum in the FOV
        """
        if (ncube is None):
            ncube = self.cube

        self.catch_line = ncube.mean(axis=(1,2))
        return self.catch_line

    def baseline(self, line=None, deduct=True, win=None, unit="km/s", order=1, full=False):
        """
        Perform a baseline fit (currently use np.polyfit)
        return a rms and
        a deducted line
        """
        if(win is None):
            raise ValueError("In function baseline() parameter \"win\" must be specific.")
        if(line is None):
            line = self.catch_line

        tmp_line = line.copy()
        line_velo = self._line_velo(line=line)
        tmp_line[self._set_win(line_velo, win)] = np.nan
        idx = np.isfinite(line_velo) & np.isfinite(tmp_line)
        parr = np.polyfit(line_velo[idx], tmp_line[idx], order)
        rms = np.sqrt(np.mean((np.polyval(parr, line_velo)[idx]-tmp_line[idx])**2))

        if(full):
            self.catch_line = line-np.polyval(parr, line_velo)
            return rms, self.catch_line, np.polyval(parr,line_velo)
        elif(deduct):
            self.catch_line = line-np.polyval(parr, line_velo)
            return rms, self.catch_line
        else:
            return rms

    def _line_velo(self, line=None, header=None, unit="km/s"):
        """
        Get velocity distribution from header
        """
        if (line is None):
            line = self.catch_line
        if (header is None):
            header = self.header
        line_ind = np.arange(line.size)
        proj = wcs(self.header)
        line_ax = ["spectral"]
        # sub is different from the original one
        line_wcs = proj.sub(line_ax)
        line_velo = line_wcs.wcs_pix2world(line_ind, 0)
        self.catch_line_velo = line_velo[0]/1000. # convert to unit "km/s"
#        print("line_velo shape is:")
#        print(line_velo)
        return self.catch_line_velo

    def _set_win(self, velo, vrange, unit="km/s", inverse=False):
        """
        Set a window range on a spectrum
        """
        if (len(vrange)%2!=0):
            raise ValueError("Velocity range must be an even number!")
        vrange_list = list(vrange)
        win_velo = np.full(velo.shape, False)
        while(len(vrange_list)>0):
            vmin = vrange_list[0]
            vmax = vrange_list[1]
            #print(vmin, vmax)
            if (vmin>vmax):
                vtmp = vmin
                vmin = vmax
                vmax = vtmp
            #remove the first two elements.
            vrange_list.pop(0)
            vrange_list.pop(0)
            win_velo = np.logical_or(win_velo, np.logical_and(velo<=vmax, velo>=vmin))
            #print(win_velo.any())

        if (inverse):
            win_velo = np.logical_not(win_velo)

        return win_velo

    def get_grid_spec(self, xsize, ysize, weight="distance"):
        """
        To get the spectra in 2D grids
        """
        (nz, ny, nx) = self.cube.shape
        if (xsize>nx or ysize>ny):
            raise ValueError("Grid size must be smaller than the original pixel size!")
        numx = np.floor(nx/xsize)
        numy = np.floor(ny/ysize)
        self.grid_spec = np.zeros((nz, ysize, xsize))
        self.grid_header = self.header.copy()
        if (self.grid_header["NAXIS"]==4):
            try:
                new_header.__delitem__("NAXIS4")
                new_header.__delitem__("CTYPE4")
                new_header.__delitem__("CRPIX4")
                new_header.__delitem__("CDELT4")
                new_header.__delitem__("CRVAL4")
                new_header.__delitem__("CROTA4")
            except KeyError:
                pass
            finally:
                new_header["NAXIS"] = 3


        pass


#    def channel_maps():


#=================Class Cube======================

#=================Function mytrim=================
#=================Function mytrim=================



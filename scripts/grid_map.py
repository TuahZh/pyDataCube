#!/usr/bin/env python
#filename=grid_map
#---------------------------------------by Gao John
#This program is designed to generate a grid of 12CO spectra 
#and overlaped with radio countour.
#The fits file name is "../kes41_12coall_MEANeff_based.fits"
#which can be changed if necessary

from __future__ import division                 #a future version of division
from pylab import *                             #contains many useful modules like numpy
import pyfits                                   #operating fits file need this module
#import pywcs                                    #it provides a method to convert coordinate
from sys import argv                            #for different files
#from string import atoi,atof                    #this can transform a string to int or float
from matplotlib.axes import Axes                #a class contains most of the figure elements
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
                                                #just as the name says
import pywcsgrid2                               #This is for displaying astronomical fits image used with matplotlib
import mpl_toolkits.axes_grid1.axes_grid as axes_grid
                                                #some useful classes and functions for grid

#function setup_axes()
def setup_axes(fig, header):
    
    ax_t = pywcsgrid2.subplot(111, header=header)
    
    #gh = pywcsgrid2.GridHelper(wcs=header)
    #gh.locator_params(nbins=4)

    g = axes_grid.ImageGrid(fig, 111, 
            nrows_ncols=(ny,nx), 
            ngrids=None, direction='row', 
            axes_pad=0.0, 
            add_all=True, 
            share_all=True, 
            aspect=True, 
            label_mode='L', 
            cbar_mode=None) 

    return g, ax_t

#function lim_up_down
def lim_up_down(vcut1,vcut2):
    if (vcut1>vcut2):
        v_temp = vcut1
        vcut1 = vcut2
        vcut2 = v_temp

    vn = ((arange(co_header['NAXIS3'])-co_header['CRPIX3'])*co_header['CDELT3']+co_header['CRVAL3'])/1000.
    k = 0
    kcut1 = 0
    kcut2 = int(co_header['CRPIX3'])
    while (k<nn-1):
        if (vn[k]<=vcut1 and vn[k+1]>=vcut1):# or (vn[k]>=-200 and vn[k+1]<=-200)):
            kcut1 = k
            break
        k+=1
    #k = 0
    while (k<nn):
        if (vn[k]<=vcut2 and vn[k+1]>=vcut2):
            kcut2 = k
            break
        k+=1

    return kcut1, kcut2


#===============================================
#file name:
fname = argv[1]
#===============================================
#================================================
#co file: ctb37b
file_co = pyfits.open(fname)

co_data = file_co[0].data
co_header = file_co[0].header
(nn, ny, nx) = co_data.shape
#===============================================
#=================================================
#radio file: ctb37b
file_radio = pyfits.open('..//most//MOS336p5_smaller_fk5.fits')
radio_data = file_radio[0].data
radio_header = file_radio[0].header

#==================================================

#test
#nx_test = 3
#ny_test = 4
#nx = nx_test
#ny = ny_test
#===============================================
#the velo-window of the spectra
vcut1 = -66
vcut2 = -35
#==============================================
#draw image
images = []
fig = figure(1, figsize=(nx,ny))
print 'it is preparing axes, please wait...'
g, ax_t = setup_axes(fig, co_header)
print 'ready?'
kcut1, kcut2 = lim_up_down(vcut1,vcut2)
#print kcut1,kcut2
nk = kcut2-kcut1
ml = 0
ii = 0
jj = ny-1
print 'go!'
for i, ax in enumerate(g):
    #nomalization
    value_max = max(co_data[:,jj,ii])
    co_data[kcut1:kcut2,jj,ii] /= value_max
  
    co_spec = zeros(nk)
    co_spec = co_data[kcut1:kcut2,jj,ii]*(nk-1)#multiply the maximum value
    im = ax.plot(co_spec)
    print i,'(',ii,',',jj,')'#if it didn't react anymore I will feel booring
    ax.set_ylim((ml,ml+nk-1))
    ax.set_xlim((0,nk-1))
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.set_aspect('equal')
    #ax.set_alpha(0.5)
    #ax.text(nk/2,0,'(%d,%d)'%(ii,jj))

    images.append(im)
    #change ii, jj
    if (ii>=nx-1):
        jj -= 1
        ii = 0
    else:
        ii += 1

#ax_t = pywcsgrid2.subplot(111, header=co_header)
#ax_t.set_alpha(0.5)
#draw contour
cont = ax_t[radio_header].contour(radio_data, levels=arange(0.,1.2,0.1), colors="y")
ax_t.set(xlim=(0,nx-1), ylim=(0,ny-1))
ax_t.grid(True)
ax_t.set_frame_on(False)
#ax_t.set_xlabel("Right Ascension (J2000)")
#ax_t.set_ylabel("Declination (J2000)")
#ax_t.axvline(x=249.8)
#show()
#ax_t.axis[:].toggle(all=False)
savefig("grid_12co_spec.eps", transparent=True, bbox_inches="tight")


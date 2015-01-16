#!/usr/bin/env python
#filename=area_12_x
#---------------------------------------by Gao John
#this file is to generate the integretion of the cube in CO observation from velocity argv[1] to argv[2]

from __future__ import division
from pylab import *
import pyfits
from sys import argv
from string import atoi,atof
from matplotlib.axes import Axes#
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable#
import pywcsgrid2

#function setup_axes()
def setup_axes():
    ax = pywcsgrid2.subplot(111, header=co_header)

    #add colorbar axes
    #divider = make_axes_locatable(ax)
    #cax = divider.new_horizontal("5%", pad=0.1, axes_class=Axes)
    #fig.add_axes(cax)

    return ax#, cax

#================================================
#co file: kes 41
file_co = pyfits.open('..//kes41_12coall_MEANeff_based.fits')

co_data = file_co[0].data
co_header = file_co[0].header
(nn, ny, nx) = co_data.shape

#=================================================
#radio file: kes 41
file_radio = pyfits.open('..//most//MOS336p5_smaller_fk5.fits')
radio_data = file_radio[0].data
radio_header = file_radio[0].header

#==================================================
#x-ray file: kes 41
file_x = pyfits.open('..//xmm//adapt-2000-6000_smaller.fits')
x_data = file_x[0].data
x_header = file_x[0].header

#=================================================
#grid helper
#grid_helper = pywcsgrid2.GridHelper(wcs=co_header)

#AxesGrid to display
#fig = figure(1)
#grid = AxesGrid(fig,)
vcut1 = atof(argv[1])
vcut2 = atof(argv[2])
if (vcut1>vcut2):
    vcut1 = atof(argv[2])
    vcut2 = atof(argv[1])

vn = ((arange(nn)-co_header['CRPIX3'])*co_header['CDELT3']+co_header['CRVAL3'])/1000.
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

#ra_x = (arange(nx)-co_header['CRPIX1'])*co_header['CDELT1']+co_header['CRVAL1']
#dec_y = (arange(ny)-co_header['CRPIX2'])*co_header['CDELT2']+co_header['CRVAL2']
co_after_cut = zeros((ny,nx))
for i in range(nx):
    for j in range(ny):
        #ns = 0.
        for k in range(nn)[kcut1:kcut2+1]:
            if (isnan(co_data[k,j,i])):
                continue
            co_after_cut[j,i] += co_data[k,j,i]
        #    ns += 1.
        #if (ns == 0.):
        #    co_after_cut[i,j] = nan
        #co_after_cut[i,j] /= ns

#print ra_x[0],ra_x[nx-1],dec_y[0],dec_y[ny-1]
#axis([ra_x[0],ra_x[nx-1],dec_y[0],dec_y[ny-1]])
#imshow(co_after_cut)
#ax = pywcsgrid2.subplot(111, header=co_header)

#============================================
#prepare figure & axes
fig = figure(1)
ax = setup_axes()

#draw image

im_x = ax[x_header].imshow(x_data, cmap=cm.hot_r, origin='lower')
#im.set_clim()

#draw contour
cont = ax[radio_header].contour(radio_data, levels=arange(0.,1.2,0.2), colors="w", alpha=0.5)
im = ax.contour(co_after_cut, levels=arange(700.,560.,-50.), colors='g', alpha=0.2)
for col in cont.collections:
    col.set_linewidth(0.5)

#cbar = colorbar(im, cax=cax)

#adjust colorbar ticks and add levels for contour lines
#cbar.set_ticks()
#cbar.add_lines(cont)
#labels
#cax.set_ylabel("Jy/Beam")

ax.set(xlim=(0,nx-1), ylim=(0,ny-1))
ax.set_xlabel("Right Ascension (J2000)")
ax.set_ylabel("Declination (J2000)")
#show()

#imshow(co_after_cut)
#colorbar()
#cont = ax[radio_header].contour(radio_data)
fnames = 'area_'+'%d'%vcut1+'_%d'%vcut2+'_12_x_v1.eps'
#ax.set_rasterization_zorder(2.1)
#cax.set_rasterization_zorder(2.1)
savefig(fnames, bbox_inches="tight")

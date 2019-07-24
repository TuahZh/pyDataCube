#!/usr/bin/env python
#filename=area_x
#---------------------------------------by Gao John

from __future__ import division
from pylab import *
import pyfits
#from sys import argv
#from string import atoi,atof
from matplotlib.axes import Axes#
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable#
import pywcsgrid2

#function setup_axes()
def setup_axes():
    ax = pywcsgrid2.subplot(111, header=x_header)

    #add colorbar axes
    divider = make_axes_locatable(ax)
    cax = divider.new_horizontal("5%", pad=0.1, axes_class=Axes)
    fig.add_axes(cax)

    return ax, cax

#================================================
#x file: kes 41
file_x = pyfits.open('..//xmm//adapt-2000-6000_smaller.fits')

x_data = file_x[0].data
x_header = file_x[0].header
(ny,nx) = x_data.shape
#=================================================
#radio file: kes 41
file_radio = pyfits.open('..//most//MOS336p5_smaller_fk5.fits')
radio_data = file_radio[0].data
radio_header = file_radio[0].header

#==================================================
#============================================
#prepare figure & axes
fig = figure(1)
ax, cax = setup_axes()

#draw image
#x_data=x_data[5:ny-1,0:nx-1]
im = ax.imshow(x_data, cmap=cm.Purples_r, origin="lower", interpolation='none')
#im.set_clim()

#draw contour
cont = ax[radio_header].contour(radio_data, levels=arange(0.,1.2,0.1), colors="y", alpha=0.5)

for col in cont.collections:
    col.set_linewidth(0.5)

cbar = colorbar(im, cax=cax)

#adjust colorbar ticks and add levels for contour lines
#cbar.set_ticks()
#cbar.add_lines(cont)
#labels
cax.set_ylabel("Jy/Beam")

#ax.set(xlim=(0,nx-1), ylim=(5,ny-1))
ax.set_xlabel("Right Ascension (J2000)")
ax.set_ylabel("Declination (J2000)")
#show()

#imshow(co_after_cut)
#colorbar()
#cont = ax[radio_header].contour(radio_data)
fnames = 'area_x.eps'
#ax.set_rasterization_zorder(2.1)
#cax.set_rasterization_zorder(2.1)
savefig(fnames, bbox_inches="tight")

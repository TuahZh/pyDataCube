#!/usr/bin/env python
#filename=gaussian_masking
#---------------------------------------by Gao John

from __future__ import division
from pylab import *
import pyfits
from scipy import optimize
import pywcsgrid2
from matplotlib.axes import Axes#
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable#

def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*exp(
                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    total = data.sum()
    X, Y = indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    col = data[:, int(y)]
    width_x = sqrt(abs((arange(col.size)-y)**2*col).sum()/col.sum())
    row = data[int(x), :]
    width_y = sqrt(abs((arange(row.size)-x)**2*row).sum()/row.sum())
    height = data.max()
    return height, x, y, width_x, width_y

def fitgaussian(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)
    errorfunction = lambda p: ravel(gaussian(*p)(*indices(data.shape)) -
                                 data)
    p, success = optimize.leastsq(errorfunction, params)
    return p

#function lim_up_down
def lim_up_down(vcut1,vcut2):
    if (vcut1>vcut2):
        v_temp = vcut1
        vcut1 = vcut2
        vcut2 = v_temp

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

    return kcut1, kcut2

#function setup_axes()
def setup_axes():
    ax = pywcsgrid2.subplot(111, header=co_header)

    #add colorbar axes
    divider = make_axes_locatable(ax)
    cax = divider.new_horizontal("5%", pad=0.1, axes_class=Axes)
    fig.add_axes(cax)

    return ax, cax


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

#vcut
vcut1 = -65
vcut2 = -50
kcut1, kcut2 = lim_up_down(vcut1,vcut2)
co_after_cut = zeros((ny,nx))
for i in range(nx):
    for j in range(ny):
        #ns = 0.
        for k in range(nn)[kcut1:kcut2+1]:
            if (isnan(co_data[k,j,i])):
                continue
            co_after_cut[j,i] += co_data[k,j,i]


params = fitgaussian(co_after_cut)
fit = gaussian(*params)

#============================================
#prepare figure & axes
fig = figure(1)
ax, cax = setup_axes()

#draw image
im = ax.imshow(fit(*indices(co_after_cut.shape)), cmap=cm.ocean_r, origin="lower", interpolation='none')
#im.set_clim()

cbar = colorbar(im, cax=cax)

#adjust colorbar ticks and add levels for contour lines
#cbar.set_ticks()
#cbar.add_lines(cont)
#labels
cax.set_ylabel("Jy/Beam")

ax.set(xlim=(0,nx-1), ylim=(0,ny-1))
ax.set_xlabel("Right Ascension (J2000)")
ax.set_ylabel("Declination (J2000)")
#show()

#imshow(co_after_cut)
#colorbar()
#cont = ax[radio_header].contour(radio_data)
fnames = 'gaussian_masking.eps'
#ax.set_rasterization_zorder(2.1)
#cax.set_rasterization_zorder(2.1)
savefig(fnames, bbox_inches="tight")
print params

#============================================
#prepare figure & axes
fig = figure(2)
ax, cax = setup_axes()

#draw image
im = ax.imshow(co_after_cut-fit(*indices(co_after_cut.shape)), cmap=cm.ocean_r, origin="lower", interpolation='none')
#im.set_clim()

#draw contour
cont = ax[radio_header].contour(radio_data, levels=arange(0.,1.2,0.1), colors="w", alpha=0.5)

for col in cont.collections:
    col.set_linewidth(0.5)

cbar = colorbar(im, cax=cax)

#adjust colorbar ticks and add levels for contour lines
#cbar.set_ticks()
#cbar.add_lines(cont)
#labels
cax.set_ylabel("Jy/Beam")

ax.set(xlim=(0,nx-1), ylim=(0,ny-1))
ax.set_xlabel("Right Ascension (J2000)")
ax.set_ylabel("Declination (J2000)")
#show()

#imshow(co_after_cut)
#colorbar()
#cont = ax[radio_header].contour(radio_data)
fnames = 'area_'+'%d'%vcut1+'_%d'%vcut2+'_12_gmasking.eps'
#ax.set_rasterization_zorder(2.1)
#cax.set_rasterization_zorder(2.1)
savefig(fnames, bbox_inches="tight")


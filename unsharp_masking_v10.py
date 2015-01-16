#!/usr/bin/env python
#filename=unsharp_masking_v10
#---------------------------------------by Gao John

from __future__ import division
from pylab import *
import pyfits
import pywcsgrid2
import scipy.ndimage as ni
from matplotlib.axes import Axes#
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable#

def load_co():
    f = pyfits.open("..//kes41_12coall_MEANeff_based.fits")
    return f[0].header,f[0].data

def load_radio():
    f = pyfits.open("..//most//MOS336p5_smaller_fk5.fits")
    return f[0].header,f[0].data

def smooth_data(data):
    smoothed = ni.gaussian_filter(data.astype('d'),36)
    return smoothed

def setup_axes(fig, header):
    ax = pywcsgrid2.subplot(111, header=header)

    #add colorbar axes
    divider = make_axes_locatable(ax)
    cax = divider.new_horizontal("5%", pad=0.1, axes_class=Axes)
    fig.add_axes(cax)

    return ax, cax

#function lim_up_down
def lim_up_down(vcut1,vcut2):
    if (vcut1>vcut2):
        v_temp = vcut1
        vcut1 = vcut2
        vcut2 = v_temp

    #nn = co_header['NAXIS3']

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

def area_v(vcut1,vcut2):
    kcut1,kcut2 = lim_up_down(vcut1,vcut2)
    co_after_cut = zeros((ny,nx))
    for i in range(nx):
        for j in range(ny):
            #ns = 0.
            for k in range(nn)[kcut1:kcut2+1]:
                if (isnan(co_data[k,j,i])):
                    continue
                co_after_cut[j,i] += co_data[k,j,i]
    
    return co_after_cut

co_header, co_data = load_co()
radio_header, radio_data = load_radio()
(nn,ny,nx) = co_data.shape

#========================================================
#vcut
vcut1 = -53
vcut2 = -52

data_area = area_v(vcut1,vcut2)
smoothed_data = smooth_data(data_area)
res_data = data_area-smoothed_data

fig1 = figure(1)
ax1, cax1 = setup_axes(fig1, co_header)
im1 = ax1.imshow(smoothed_data, cmap=cm.ocean_r, origin='lower', interpolation='none')
cbar1 = colorbar(im1, cax=cax1)

axis = ax1["gal"].new_floating_axis(1, 0.)
axis.toggle(all=False, label=True)
axis.label.set_text(r"$\hspace{4}b=0^{\circ}$")
ax1.axis["b=0"] = axis

fname1 = 'smoothed12_v10.eps'
savefig(fname1, bbox_inches='tight')


fig2 = figure(2)
ax2, cax2 = setup_axes(fig2, co_header)
im2 = ax2.imshow(res_data, cmap=cm.ocean_r, origin='lower', interpolation='none')
im2.set_clim(vmin=0.)
cont2 = ax2[radio_header].contour(radio_data, 
        levels=arange(0.,1.2,0.1), colors='r')
cbar2 = colorbar(im2, cax=cax2)
fname2 = 'res12_data_%d_%d.eps'%(vcut1,vcut2)
savefig(fname2, bbox_inches='tight')


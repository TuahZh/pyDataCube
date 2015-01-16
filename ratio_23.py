#!/usr/bin/env python
#filename=ratio_23
#---------------------------------------by Gao John

from __future__ import division
from pylab import *
import pyfits
import pywcsgrid2
from matplotlib.axes import Axes#
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable#

def load_12co():
    f = pyfits.open("..//kes41_12coall_MEANeff_based.fits")
    return f[0].header,f[0].data

def load_radio():
    f = pyfits.open("..//most//MOS336p5_smaller_fk5.fits")
    return f[0].header,f[0].data

def load_13co():
    f = pyfits.open("..//kes41_13coall_MEANeff_based.fits")
    return f[0].header,f[0].data

def setup_axes(fig, header):
    ax = pywcsgrid2.subplot(111, header=header)

    #add colorbar axes
    divider = make_axes_locatable(ax)
    cax = divider.new_horizontal("5%", pad=0.1, axes_class=Axes)
    fig.add_axes(cax)

    return ax, cax

#function lim_up_down
def lim_up_down(vcut1,vcut2,co_header):
    if (vcut1>vcut2):
        v_temp = vcut1
        vcut1 = vcut2
        vcut2 = v_temp

    nn = co_header['NAXIS3']

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

def area_v(vcut1,vcut2,co_data,co_header,nn):
    kcut1,kcut2 = lim_up_down(vcut1,vcut2,co_header)
    co_after_cut = zeros((ny,nx))
    for i in range(nx):
        for j in range(ny):
            #ns = 0.
            for k in range(nn)[kcut1:kcut2+1]:
                if (isnan(co_data[k,j,i])):
                    continue
                co_after_cut[j,i] += co_data[k,j,i]
    
    return co_after_cut

def ratio_get(data12,data13):
    #rms
    rms12_s = 0.2*6
    rms13_s = 0.1*6
    kcut1,kcut2 = lim_up_down(vcut1,vcut2,co12_header)
    rms12 = rms12_s*abs(kcut2-kcut1)
    kcut1,kcut2 = lim_up_down(vcut1,vcut2,co13_header)
    rms13 = rms13_s*abs(kcut2-kcut1)

    ratio_data = zeros((ny,nx))
    for i in range(nx):
        for j in range(ny):
            if (abs(data12[j,i])<rms12 or abs(data13[j,i])<rms13):
                ratio_data[j,i] = 0
                continue
            ratio_data[j,i] = data12[j,i]/data13[j,i]

    return ratio_data

co12_header, co12_data = load_12co()
co13_header, co13_data = load_13co()
radio_header, radio_data = load_radio()
(nn,ny,nx) = co12_data.shape
(nn_p,ny_p,nx_p) = co13_data.shape
if (ny!=ny_p or nx!=nx_p):
    print "Please make sure 13 and 12 have the same size of their datas!"
    exit(1)


#========================================================
#vcut
vcut1 = -55
vcut2 = -50

data_area12 = area_v(vcut1,vcut2,co12_data,co12_header,nn)
data_area13 = area_v(vcut1,vcut2,co13_data,co13_header,nn_p)
ratio_data = ratio_get(data_area12,data_area13)

fig1 = figure(1)
ax1, cax1 = setup_axes(fig1, co12_header)
im1 = ax1.imshow(ratio_data, cmap=cm.ocean_r, origin='lower', interpolation='none')
im1.set_clim(vmin=2.,vmax=5.)
cont1 = ax1[radio_header].contour(radio_data, levels=arange(0.,1.2,0.1), colors='w')
cbar1 = colorbar(im1, cax=cax1)
fname1 = 'ratio_%d_%d.eps'%(vcut1,vcut2)
savefig(fname1, bbox_inches='tight')


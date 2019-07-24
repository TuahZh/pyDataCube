#!/usr/bin/env python
#filename=P-V_12
#---------------------------------------by Gao John

from __future__ import division
from pylab import *
import pyfits
import pywcsgrid2
import mpl_toolkits.axisartist as axisartist
from mpl_toolkits.axes_grid import make_axes_locatable
import matplotlib.colors as mcolors
import pywcs

#function lim_up_down
def lim_up_down(vcut1,vcut2):
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

def co_area(kcut1,kcut2):
    co_after_cut = zeros((ny,nx))
    for i in range(nx):
        for j in range(ny):
            #ns = 0.
            for k in range(nn)[kcut1:kcut2+1]:
                if (isnan(co_data[k,j,i])):
                    continue
                co_after_cut[j,i] += co_data[k,j,i]
    
    return co_after_cut

def setup_axes(fig, header):
    
    ax0 = pywcsgrid2.subplot(111, wcs=header)
    divider = make_axes_locatable(ax0)

    gh1 = pywcsgrid2.GridHelperSimple(wcs=header, axis_nums=[0, 2])
    ax_v = divider.new_vertical(4, pad=0.1, sharex=ax0,
                                axes_class=pywcsgrid2.Axes,
                                grid_helper=gh1)
    fig.add_axes(ax_v)

    gh2 = pywcsgrid2.GridHelperSimple(wcs=header, axis_nums=[2, 1])
    ax_h = divider.new_horizontal(4, pad=0.1, sharey=ax0,
                                axes_class=pywcsgrid2.Axes,
                                grid_helper=gh2)

    fig.add_axes(ax_h)


    ax_h.axis["left"].toggle(label=False, ticklabels=False)
    ax_v.axis["bottom"].toggle(label=False, ticklabels=False)


    return ax0, ax_v, ax_h

#co_file===================================================
co_file = pyfits.open("..//kes41_12coall_MEANeff_based.fits")
co_header = co_file[0].header
co_data = co_file[0].data
(nn, ny, nx) = co_data.shape
co_wcs = pywcs.WCS(co_header)

#radio_file===============================================
file_radio = pyfits.open('..//most//MOS336p5_smaller_fk5.fits')
radio_data = file_radio[0].data
radio_header = file_radio[0].header

#vcut
vcut1=-55
vcut2=-45

kcut1, kcut2 = lim_up_down(vcut1,vcut2)

#pick x, y=========================================
pick_x = 249.7
pick_y = -46.97
co_wcs_xy = co_wcs.sub([1,2])
[[pick_i,pick_j]] = co_wcs_xy.wcs_sky2pix([[pick_x,pick_y]], 0)
print pick_i, pick_j

norm = mcolors.Normalize()
cmap = cm.ocean_r

#co_file[0].data = co_file[0].data**2

vv0  = co_area(kcut1,kcut2)
vv1  = co_data[:,pick_j,:]
vv2  = co_data[:,:,pick_i].transpose()


fig = figure(1, figsize=(12, 12))
ax0, ax_v, ax_h = setup_axes(fig, co_header)

im0 = ax0.imshow(vv0, origin="lower", 
                 cmap=cmap)

cont0 = ax0[radio_header].contour(radio_data,levels=arange(0.,1.2,0.1),
        colors="w", alpha=0.5)

plt0 = ax0["fk5"].plot([pick_x], [pick_y], "w+", ms=9, mew=2, zorder=3)

im_v = ax_v.imshow(vv1, origin="lower", aspect="auto",
            cmap=cmap)

cont_v = ax_v.contour(vv1, colors='w', levels=arange(1.,6.,1.), alpha=0.5)

im_h = ax_h.imshow(vv2, origin="lower", aspect="auto",
            cmap=cmap)

cont_h = ax_h.contour(vv2, levels=arange(1.,6.,1.), colors='w', alpha=0.5)

for col in cont_v.collections:
    col.set_linewidth(0.5)
for col in cont_h.collections:
    col.set_linewidth(0.5)

ax0.grid(True)
ax_h.grid(True)
ax_v.grid(True)

#norm.vmin=0.1
#norm.vmax=5**2
#for ax in [ax0, ax_h, ax_v]:
#    ax.images[0].changed()

gh1 = ax_v.get_grid_helper()
gh1.set_ticklabel2_type("absval", scale=0.001, nbins=5)

gh2 = ax_h.get_grid_helper()
gh2.set_ticklabel1_type("absval", scale=0.001, nbins=5)

ax_h.axis["bottom"].label.set_text(r"$v_{\mathrm{LSR}}$ [km s$^{-1}$]")
ax_v.axis["left"].label.set_text(r"$v_{\mathrm{LSR}}$ [km s$^{-1}$]")

ax_h.set(xlim=(kcut1,kcut2+1))
ax_v.set(ylim=(kcut1,kcut2+1))

fname_o = "P-V_area.eps"
savefig(fname_o, bbox_inches='tight')

#
if 0:
    gh1 = pywcsgrid2.GridHelperSimple(wcs=co_header, axis_nums=[2, 0])
    gh2 = pywcsgrid2.GridHelperSimple(wcs=co_header, axis_nums=[2, 1])
    fig_n = figure(2)
    ax1 = pywcsgrid2.subplot(211, grid_helper=gh1)
    im1 = ax1.imshow(vv1.transpose(), origin='lower', aspect='auto')
    ax2 = pywcsgrid2.subplot(212, grid_helper=gh2)
    im2 = ax2.imshow(vv2, origin='lower', aspect='auto')
    fname_o_new = 'P-V_%f_%f.eps'%(pick_x,pick_y)
    ax1.set(xlim=(kcut1,kcut2+1))
    ax2.set(xlim=(kcut1,kcut2+1))
    cont1 = ax1.contour(vv1.transpose(), levels=arange(2.,6.,2.), colors='w',
            linewidths=0.5)
    cont2 = ax2.contour(vv2, levels=arange(2.,6.,2.), colors='w', 
            linewidths=0.5)
    savefig(fname_o_new, bbox_inches='tight')

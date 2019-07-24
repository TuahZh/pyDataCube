#!/usr/bin/env python
#filename=channel_maps18_v1
#---------------------------------------by Gao John
#this file is to generate the integretion of the cube in CO observation from velocity argv[1] to some velocity by 1km/s

from __future__ import division
from pylab import *
import pyfits
#import pywcs
from sys import argv
from string import atoi,atof
from matplotlib.axes import Axes#
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable#
import pywcsgrid2
import mpl_toolkits.axes_grid1.axes_grid as axes_grid
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.colors as mcolors
use_path_effect = True
try:
    from matplotlib.patheffects import withStroke
except ImportError:
    use_path_effect = False

#define class Velo
#class Velo(object):
#    def __init__(self, header):
#        wcs = pywcs.WCS(header)
#        self.wcs_vel = wcs.sub([3])
#
#    def to_vel(self, p):
#        v = self.wcs_vel.wcs_pix2sky([[p]], 0)
#        return v[0][0]

#function setup_axes()
def setup_axes(fig, header):
    
    gh = pywcsgrid2.GridHelper(wcs=header)
    gh.locator_params(nbins=3)

    g = axes_grid.ImageGrid(fig, 111, 
            nrows_ncols=(2,3), 
            ngrids=None, direction='row', 
            axes_pad=0.02, 
            add_all=True, 
            share_all=True, 
            aspect=True, 
            label_mode='L', 
            cbar_mode=None, 
            axes_class=(pywcsgrid2.Axes, dict(grid_helper=gh)))

    #make colorbar
    ax = g[-1]
    cax = inset_axes(ax, 
            width="8%", 
            height="100%", 
            loc=3,
            bbox_to_anchor=(1.02,0,1,1),
            bbox_transform=ax.transAxes,
            borderpad=0.
            )


    return g, cax

#================================================
#co file: ctb37b
file_co = pyfits.open('..//kes41_c18oall_MEANeff_based.fits')

co_data = file_co[0].data
co_header = file_co[0].header
(nn, ny, nx) = co_data.shape

#=================================================
#radio file: ctb37b
file_radio = pyfits.open('..//most//MOS336p5_smaller_fk5.fits')
radio_data = file_radio[0].data
radio_header = file_radio[0].header

#==================================================

#grid helper
#grid_helper = pywcsgrid2.GridHelper(wcs=co_header)

#AxesGrid to display
#fig = figure(1)
#grid = AxesGrid(fig,)

vcut1 = atof(argv[1])
vcut_for_name = vcut1
#vcut2 = atof(argv[2])
#if (vcut1>vcut2):
#    vcut1 = atof(argv[2])
#    vcut2 = atof(argv[1])
#
vn = ((arange(nn)-co_header['CRPIX3'])*co_header['CDELT3']+co_header['CRVAL3'])/1000.


#ra_x = (arange(nx)-co_header['CRPIX1'])*co_header['CDELT1']+co_header['CRVAL1']
#dec_y = (arange(ny)-co_header['CRPIX2'])*co_header['CDELT2']+co_header['CRVAL2']

#velo
#vel = Velo(co_header)

#co_after_cut = zeros((nx,ny))
#for i in range(nx):
#    for j in range(ny):
#        #ns = 0.
#        for k in range(nn)[kcut1:kcut2+1]:
#            if (isnan(co_data[k,j,i])):
#                continue
#            co_after_cut[i,j] += co_data[k,j,i]
#        #    ns += 1.
        #if (ns == 0.):
        #    co_after_cut[i,j] = nan
        #co_after_cut[i,j] /= ns

#print ra_x[0],ra_x[nx-1],dec_y[0],dec_y[ny-1]
#axis([ra_x[0],ra_x[nx-1],dec_y[0],dec_y[ny-1]])
#imshow(co_after_cut)
#ax = pywcsgrid2.subplot(111, header=co_header)

#============================================
#prepare figure & axes
fig = figure(1, figsize=(12,16), dpi=70)
g, cax = setup_axes(fig, co_header)

#draw image
i = 0
#dxy = kcut1
#nxy = 3*3
#cmap = cm.ocean_r
#cmap = cm.jet
cmap = cm.Greys
#norm =  mcolors.Normalize()
images = []
#start_channel = i*nxy+dxy

del_v = 2

for i, ax in enumerate(g):
   # channel_number = start_channel + i
   # channel = co_data[channel_number]
    k = 0
    kcut1 = 0
    kcut2 = int(co_header['CRPIX3'])
    while (k<nn-1):
        if (vn[k]<=vcut1 and vn[k+1]>=vcut1):# or (vn[k]>=-200 and vn[k+1]<=-200)):
            kcut1 = k
            break
        k+=1
    k = 0
    while (k<nn):
        if (vn[k]<=vcut1+del_v and vn[k+1]>=vcut1+del_v):
            kcut2 = k
            break
        k+=1

    co_per1 = zeros((ny,nx))
    for ii in range(nx):
        for jj in range(ny):
#        #ns = 0.
            for kk in range(nn)[kcut1:kcut2+1]:
                if (isnan(co_data[kk,jj,ii])):
                    continue
                co_per1[jj,ii] += co_data[kk,jj,ii]

    im = ax.imshow(co_per1, origin="lower",  cmap=cmap, vmin=6.0, interpolation="none") #3sigma
    #draw contour
    cont = ax[radio_header].contour(radio_data, linewidth=0.5, levels=arange(0.,1.2,0.1),
            colors="g", alpha=0.2)
    v = vcut1
    #maser point
    m_alpha = 15.*(16.+(38.+52.15/60.)/60)
    m_delta = -(46.+(56.+16.1/60)/60)
    plt = ax["fk5"].plot([m_alpha], [m_delta], "m+", ms=9, mew=2, zorder=3)

    t = ax.add_inner_title(r"$%4.1f$" % (v), loc=3, frameon=False)
    if use_path_effect:
        t.txt._text.set_path_effects([withStroke(foreground='w',
            linewidth=1)])

    images.append(im)
    vcut1 += del_v

#label with velocities
#for i, ax in enumerate(g):
#    channel_number = start_channel+i
#    v = vn[channel_number]
#    t = ax.add_inner_title(r"$v=%4.1f$ km s$^{-1}$" % (v), loc=2, frameon=False)
#    if use_path_effect:
#        t.txt._text.set_path_effects([withStroke(foreground='w',
#            linewidth=3)])

#make colorbar
cb =  colorbar(im, cax=cax)
cb.set_label("Jy/Beam")
#cb.set_ticks([0,1,2,3])

#adjust norm
#norm.vmin = -0.1
#norm.vmax = 5.0
#for im in images:
#    im.changed()

savefig("co18_channel_maps_%d"%vcut_for_name+"_v1.eps", bbox_inches="tight")
#
#im = ax.imshow(co_after_cut, cmap=cm.ocean_r, origin="lower", interpolation="nearest")
#im.set_clim()

#draw contour
#cont = ax[radio_header].contour(radio_data, levels=arange(0.,3.2,0.2), colors="w", alpha=0.5)

#for col in cont.collections:
#    col.set_linewidth(0.5)

#cbar = colorbar(im, cax=cax)

#adjust colorbar ticks and add levels for contour lines
#cbar.set_ticks()
#cbar.add_lines(cont)
#labels
#cax.set_ylabel("Jy/Beam")

#ax.set(xlim=(0,38), ylim=(0,36))
#ax.set_xlabel("Right Ascension (J2000)")
#ax.set_ylabel("Declination (J2000)")
#show()

#imshow(co_after_cut)
#colorbar()
#cont = ax[radio_header].contour(radio_data)
#fnames = 'area_'+'%d'%vcut1+'_%d'%vcut2+'_12_v1.ps'
#ax.set_rasterization_zorder(2.1)
#cax.set_rasterization_zorder(2.1)
#savefig(fnames, format='ps')

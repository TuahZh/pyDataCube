#!/usr/bin/env python
#filename=channel_maps_v2
#---------------------------------------by Gao John
#this file is to generate the integretion of the cube in CO observation from velocity argv[1] to some velocity by 1km/s

from __future__ import division
from pylab import *
import pyfits
#import pywcs
from kapteyn import wcs
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


import gaoyuan

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
    gh.locator_params(nbins=2)

    g = axes_grid.ImageGrid(fig, 111, 
            nrows_ncols=(5,5), 
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
#h1 file: kes41
file_h1 = pyfits.open('../kes41_12coall_MEANeff_based.fits')

h1_data = file_h1[0].data
h1_header = file_h1[0].header
h1_cube = gaoyuan.Cube(h1_data, h1_header).cube
(nn, ny, nx) = h1_cube.shape

#=================================================
#radio file: kes41
file_radio = pyfits.open('..//most//MOS336p5_smaller.fits')
radio_data = file_radio[0].data
radio_header = file_radio[0].header

#==================================================

vcut1 = atof(argv[1])
vcut_for_name = vcut1
#vcut2 = atof(argv[2])
#if (vcut1>vcut2):
#    vcut1 = atof(argv[2])
#    vcut2 = atof(argv[1])
#
#vn = ((arange(nn)-co_header['CRPIX3'])*co_header['CDELT3']+co_header['CRVAL3'])/1000.


#============================================
#prepare figure & axes
fig = figure(1, figsize=(12,16))
g, cax = setup_axes(fig, h1_header)

#draw image
i = 0

cdict = {'red':   [(0.0,  1.0, 1.0),
                   (0.5,  0.75, 0.75),
                   (0.71,  0.5,  0.5),
#                   (0.87,  0.25, 0.25),
                   (1.0,  0.0, 0.0)],

         'green': [(0.0,  1.0, 1.0),
                   (0.5,  0.75, 0.75),
                   (0.71,  0.5,  0.5),
#                   (0.87,  0.25, 0.25),
                   (1.0,  0.0, 0.0)],

         'blue':  [(0.0,  1.0, 1.0),
                   (0.5,  0.75, 0.75),
                   (0.71,  0.5,  0.5),
#                   (0.87,  0.25, 0.25),
                   (1.0,  0.0, 0.0)]}

my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)
#dxy = kcut1
#nxy = 3*3
#cmap = cm.ocean_r
#cmap = cm.gist_ncar
#cmap = cm.Greys
#cmap = cm.jet
#cmap = cm.gist_yarg
#cmap = cm.binary
#cmap = cm.YlOrRd
cmap = my_cmap
#norm =  mcolors.Normalize()
images = []
#start_channel = i*nxy+dxy
proj = wcs.Projection(h1_header)
line_ax = proj.specaxnum
line = proj.sub(line_ax)

del_v = 1
#unit change
vcut1 = vcut1*1e3
del_v = del_v*1e3
for i, ax in enumerate(g):
   # channel_number = start_channel + i
   # channel = co_data[channel_number]
#    k = 0
#    kcut1 = 0
#    kcut2 = int(co_header['CRPIX3'])
#    while (k<nn-1):
#        if (vn[k]<=vcut1 and vn[k+1]>=vcut1):# or (vn[k]>=-200 and vn[k+1]<=-200)):
#            kcut1 = k
#            break
#        k+=1
#    k = 0
#    while (k<nn):
#        if (vn[k]<=vcut1+del_v and vn[k+1]>=vcut1+del_v):
#            kcut2 = k
#            break
#        k+=1
#
    kcut1 = int(line.topixel1d(vcut1))
    vcut2 = vcut1+del_v
    kcut2 = int(line.topixel1d(vcut2))

    co_per1 = zeros((ny,nx))
    for ii in range(nx):
        for jj in range(ny):
#        #ns = 0.
            for kk in range(nn)[kcut1:kcut2+1]:
                if (isnan(h1_cube[kk,jj,ii])):
                    print 'nan found'
                    continue
                co_per1[jj,ii] += h1_cube[kk,jj,ii]

    im = ax.imshow(co_per1, origin="lower",  cmap=cmap, 
            vmin=14, 
            vmax=56.5, 
#            vmax=co_per1.max(),
            interpolation="none") 
    #draw contour
#    cont = ax[radio_header].contour(radio_data, linewidth=0.5,
#            levels=arange(0.07,0.88,0.06), colors="c", alpha=0.8)
    cont = ax[radio_header].contour(radio_data, linewidth=0.5,
           levels=arange(0.07,0.88,0.12), colors="c", alpha=0.8)
    v = vcut1/1e3
    t = ax.add_inner_title(r"$%4.1f$" % (v), loc=3, frameon=False,
            prop={'size':16})
#maser point
    m_alpha = 15.*(16.+(38.+52.15/60.)/60)
    m_delta = -(46.+(56.+16.1/60)/60)
    plt = ax["fk5"].plot([m_alpha], [m_delta], color="#fa626e", marker="+", ms=20, mew=2, zorder=3)

    ledge=(249.8672754223264, -46.90550387875837)
    redge=(249.6, -47.046662665497005)
    lshell=(249.84241, -46.918825)
    rshell=(249.68954, -46.998512)
    cent = (0.5*(lshell[0]+rshell[0]),0.5*(lshell[1]+rshell[1]))
#column density 
    if (i==0):
        plt_line = ax["fk5"].plot([ledge[0],redge[0]],[ledge[1],redge[1]],'--', color='#ff6000')
        plt_cen = ax["fk5"].plot([cent[0]],[cent[1]], color='#ff6000', marker='x')

    if use_path_effect:
        t.txt._text.set_path_effects([withStroke(foreground='w',
            linewidth=1)])

    ax.axis["bottom"].major_ticklabels.set(fontsize=18)
    ax.axis["bottom"].label.set(fontsize=18)
    ax.axis["left"].major_ticklabels.set(fontsize=18)
    ax.axis["left"].label.set(fontsize=18)
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
cb.set_label("K",fontsize=18)
cb.set_ticks([15,25,35,45,55])
for t in cb.ax.get_yticklabels():
    t.set_fontsize(18)


#cb.set_ticks([28,36,44,52])

#adjust norm
#norm.vmin = 125
#norm.vmax = 300
#for im in images:
#    im.changed()

savefig("co12_channel_maps_%d"%vcut_for_name+"_v3_1.eps",
        bbox_inches="tight", dpi=300, format="eps")



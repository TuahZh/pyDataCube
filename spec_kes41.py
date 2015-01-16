#!/usr/bin/env python
#filename=spec_kes41
#---------------------------------------by Gao John

from __future__ import division
from pylab import *
import pyfits

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

#---------------------------------------
v_cut1 = -150
v_cut2 = 10
#---------------------------------------
ax = subplot(111)

#file names
for n_f in ['12co','13co','c18o']:
    fname = '..//kes41_%sall_MEANeff_based.fits' % n_f
    hdulist = pyfits.open(fname)
    (nn, ny, nx) = hdulist[0].data.shape
    imdata = hdulist[0].data
    co_header = hdulist[0].header

    k_1, k_2 = lim_up_down(v_cut1, v_cut2)

    ss = zeros(nn)
    for i in range(nn):
        nu = 0.
        for j in range(nx):
            for k in range(ny):
                if isnan(imdata[i,k,j]): 
                    continue
                nu += 1
                ss[i] += imdata[i,k,j]
        ss[i] /= nu

    header =  hdulist[0].header

    index_v = ((arange(header['NAXIS3'])-header['CRPIX3'])*header['CDELT3']+header['CRVAL3'])/1000.
    #print ss
    n_ff={"12co":r"$^{12}$CO","13co":r"$^{13}$CO","c18o":r"C$^{18}$O"}
#main beam efficiency
    eta = 0.42
    ax.plot(index_v[k_1:k_2],ss[k_1:k_2]/eta, label='%s' % n_ff[n_f])

ax.set_xlim(v_cut1,v_cut2)
ax.set_xlabel('Velocity (km/s)',fontsize=18)
ax.set_ylabel(r'$T_{R}$ (K)',fontsize=18)
ax.legend(prop={'size':18})
ax.axvline(x=-45, color='k', ls=':')
for t in ax.xaxis.get_major_ticks():
    t.label.set_fontsize(18)
 
for t in ax.yaxis.get_major_ticks():
    t.label.set_fontsize(18)

savefig('spc_t.eps', bbox_inches='tight')

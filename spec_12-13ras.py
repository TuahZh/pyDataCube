#!/usr/bin/env python
#filename=spec_12-13
#---------------------------------------by Gaoyuan

#import sometings
from __future__ import division
from pylab import *
#from sys import argv
import pyfits
#import os
from scipy import interpolate

def load_12():
    f = pyfits.open("..//kes41_12coall_MEANeff_based.fits")
    return f[0].header,f[0].data

def load_13():
    f = pyfits.open("..//kes41_13coall_MEANeff_based.fits")
    return f[0].header,f[0].data

def extract_spectrum((x,y),data):
    line = (data[:,y,x]+data[:,y+1,x+1]+data[:,y,x+1]+data[:,y+1,x])/4.
    return line

def locate_xy(alpha,delta,header):
    ii = (arange(header['NAXIS1'])-header['CRPIX1'])*header['CDELT1']+header['CRVAL1']
    jj = (arange(header['NAXIS2'])-header['CRPIX2'])*header['CDELT2']+header['CRVAL2']

    for i in range(header['NAXIS1']-1):
        if ((alpha>=ii[i] and alpha<ii[i+1]) or (alpha<=ii[i] and alpha>ii[i+1])):
            x = i
        if ((delta>=jj[i] and delta<jj[i+1]) or (delta<=jj[i] and delta>jj[i+1])):
            y = i

    print ii[x],jj[y]
    return (x,y)

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

def plot_line(line,header,n_f):
    index_v = ((arange(header['NAXIS3'])-header['CRPIX3'])*header['CDELT3']+header['CRVAL3'])/1000.
    plot(index_v,line, label=r'$^{%d}CO$'% n_f)

#============================
#define alpha, delta:
alpha = 15.*(16.+(38.+52.15/60.)/60.)
delta = -(46.+(56.+16.1/60.)/60.)
print alpha,delta
#vcut:
vcut1 = -150
vcut2 = +20

#load files:
co12_header, co12_data = load_12()
co13_header, co13_data = load_13()

#12CO:
(x,y) = locate_xy(alpha,delta,co12_header)
line12 = extract_spectrum((x,y),co12_data)
plot_line(line12,co12_header,12)
print x,y

#13CO:
(x,y) = locate_xy(alpha,delta,co13_header)
line13 = extract_spectrum((x,y),co13_data)
plot_line(line13,co13_header,13)
print x,y

xlim(vcut1,vcut2)
xlabel('Velocity(km/s)')
ylabel('Intensity(Jy/Beam)')
legend()
axvline(x=-45, color='r')
axvline(x=-52, color='y')
savefig('spc_m.eps', bbox_inches='tight')

#del & res:
max12 = line12.max()
max13 = line13.max()
rat = max12/max13

if 1:
    figure(2)
#12:
    kcut1,kcut2 = lim_up_down(vcut1,vcut2,co12_header)
    line12 = line12[kcut1:kcut2+1]
    index_v = ((arange(co12_header['NAXIS3'])-co12_header['CRPIX3'])*co12_header['CDELT3']+co12_header['CRVAL3'])/1000.
    index_v = index_v[kcut1:kcut2+1]
    f12 = interpolate.interp1d(index_v,line12)
    #13:
    line13 = line13*rat
    kcut1,kcut2 = lim_up_down(vcut1,vcut2,co13_header)
    line13 = line13[kcut1:kcut2+1]
    index_v = ((arange(co13_header['NAXIS3'])-co13_header['CRPIX3'])*co13_header['CDELT3']+co13_header['CRVAL3'])/1000.
    index_v = index_v[kcut1:kcut2+1]
    f13 = interpolate.interp1d(index_v,line13)
    index_new = arange(vcut1,vcut2,0.1)

    plot(index_new,f12(index_new)-f13(index_new),color='r')
    #plot(index_new,f12(index_new),color='g')
    #plot(index_new,f13(index_new),color='y')
    xlim(vcut1,vcut2)
    xlabel('Velocity(km/s)')
    ylabel('Intensity(Jy/Beam)')
    axvline(x=-45, color='k')
    axvline(x=-52, color='k')

    savefig('spc_r.eps', bbox_inches='tight')




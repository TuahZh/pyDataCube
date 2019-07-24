#!/usr/bin/env python
#filename=area_12_v1
#---------------------------------------by Gao John
#this file is to generate the integretion of the cube in CO observation from velocity argv[1] to argv[2]
#The rms of this file is recorded in file rms_area.out

from __future__ import division
from pylab import *
import pyfits
from sys import argv
from string import atoi,atof

#================================================
#co file: kes 41
file_co = pyfits.open('..//kes41_13coall_MEANeff_based.fits')

co_data = file_co[0].data
co_header = file_co[0].header
(nn, ny, nx) = co_data.shape

#==================================================
#file to record:
vcut1 = atof(argv[1])
vcut2 = atof(argv[2])
f_rms = open("rms_area13_%d_%d.out"%(vcut1,vcut2),'w')
f_read = open('rms_kes41_13coall_MEANeff.out')


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

co_after_cut = zeros((ny,nx))
for i in range(nx):
    for j in range(ny):
        for k in range(nn)[kcut1:kcut2+1]:
            if (isnan(co_data[k,j,i])):
                continue
            co_after_cut[j,i] += co_data[k,j,i]

fnames = 'kes41_'+'%d'%vcut1+'_%d'%vcut2+'_13_v1.fits'

file_co[0].data = co_after_cut
co_header.add_history("compress to "+fnames+" by GY")
co_header['NAXIS'] = 2
co_header.__delitem__('NAXIS3')
co_header.__delitem__('CTYPE3')
co_header.__delitem__('CRPIX3')
co_header.__delitem__('CDELT3')
co_header.__delitem__('CRVAL3')

file_co.writeto(fnames)

f_rms.write(fnames+':\n')
rms_mean = 0
count = 0
for ii in range(nx):
    for jj in range(ny):
        str_rms = f_read.readline()
        ind_rms = str_rms.find('rms=')
        str_rms = str_rms[ind_rms+4:-1]#\n need -1
        rms_total = atof(str_rms)*sqrt(kcut2-kcut1)
        f_rms.write("%d\t%d\t%f\n"%(ii,jj,rms_total))
        rms_mean += rms_total
        count += 1

f_rms.write("mean\t\t%f\n"%(rms_mean/count))
f_rms.close()
f_read.close()




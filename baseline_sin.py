#!/usr/bin/env python
#filename=baseline_sin
#---------------------------------------by Gao John
#This file is designed to generate a new fits file
#that contains data from a CO (or other) spectrum observing 
#with the spectra reduce the base line of sin or cos shape.

#You need to INPUT file name of the initial file and
#the 2 window velocities in which there are emission lines.

from __future__ import division
from pylab import *
import pyfits
from sys import argv
from string import atof, atoi
from scipy import optimize
import os
import time

#==========================================
#Defination part
#class that is used by parameters of the fit function.
class Parameter:
    def __init__(self, value):
            self.value = value

    def set(self, value):
            self.value = value

    def __call__(self):
            return self.value

def fit(function, parameters, y, x = None):
    def f(params):
        i = 0
        for p in parameters:
            p.set(params[i])
            i += 1
        return y - function(x)

    if x is None: x = arange(y.shape[0])
    p = [param() for param in parameters]
    out = optimize.leastsq(f, p, full_output=1)
    covar = out[1]
    #if(sqrt(covar[0][0])>0.5 or sqrt(covar[1][1])>0.5):
    #        print ii,jj,"covar too large!"
    #        exit()

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

#function rms()
def rms(data):
    l1 = data[:kcut1,jj,ii]
    l2 = data[kcut2:,jj,ii]
    rms1 = sum(data[:kcut1,jj,ii]**2)
    rms2 = sum(data[kcut2:,jj,ii]**2)
    rms_t = sqrt((rms1+rms2)/(l1.size+l2.size))
    return rms_t

#==============================================
#file name:
fname = argv[1]
ind = fname.index('.',-5)
fname_o = fname[:ind]+'_basedsin.fits'

#to make sure there is no such a file
os.system('rm -rf '+fname_o)

#file to record rms
f_rms = open("rms_"+fname[3:ind]+".out",'w')
rms_mean = 0
count = 0

#2 velocities as a window
vcut1 = atof(argv[2])
vcut2 = atof(argv[3])

#open the original file
hdulist = pyfits.open(fname)
co_header = hdulist[0].header
co_data = hdulist[0].data

#extract the spectra
[nn,ny,nx] = co_data.shape
#test part
#ii = 17
#jj = 14
new_data = co_data#define new_data
kcut1, kcut2 = lim_up_down(vcut1,vcut2)
#print kcut1,kcut2


#initial parameters:
a = Parameter(1.)
b = Parameter(0.02)
c = Parameter(0)

k_v = range(nn)#standard python list can add items easily
k_v_after = array(k_v[:kcut1]+k_v[kcut2:])#not add the value but get a larger list

#print a(), b()
#define the function
y = lambda x: a()*sin(b()*x+c()) #a and b are types of Parameter defined above

        
vn = ((arange(nn)-co_header['CRPIX3'])*co_header['CDELT3']+co_header['CRVAL3'])/1000.

for ii in xrange(nx):
    for jj in xrange(ny):


        #rebin size
        rebin_s = 10
        rebin_data=arange(kcut1//rebin_s+(nn-1-kcut2)//rebin_s)
        for kk in range(kcut1//rebin_s-1):
            kkk = kk*rebin_s
            rebin_data[kk] = mean(co_data[kkk:kkk+rebin_s,jj,ii])

        rebin_data[kcut1//rebin_s-1] = mean(co_data[(kcut1//rebin_s-1)*rebin_s:kcut1,jj,ii])

        for kk in range((nn-1-kcut2)//rebin_s-1):
            kkk = kk*rebin_s+kcut2
            rebin_data[kcut1//rebin_s+kk] = mean(co_data[kkk:kkk+rebin_s,jj,ii])

        rebin_data[kcut1//rebin_s+(nn-1-kcut2)//rebin_s-1] = mean(co_data[kcut2+((nn-1-kcut2)//rebin_s-1)*rebin_s:nn-1,jj,ii])

        k_v_rebin = array((arange(kcut1//rebin_s)*rebin_s).tolist() + (kcut2+rebin_s*arange((nn-1-kcut2)//rebin_s)).tolist())
        #fit!
        need_data = array(co_data[:kcut1,jj,ii].tolist()+co_data[kcut2:,jj,ii].tolist())

        #print need_data.shape,'\n',k_v_after.shape
        print ii,jj
        fit(y, [a,b,c], rebin_data, x=k_v_rebin)

        #print a(),b()
        #print y(array([2000,3000]))
        #plot the result
        if 0:
            plot(rebin_data[:], 'co')
            #plot(arange(nn), y(arange(nn)), 'r-')
            show()
            raw_input()

        k_v = array(k_v)#transform to array of numpy
        base_line = y(k_v)
        #reduce the base line
        new_data[:,jj,ii] = co_data[:,jj,ii]-base_line
        f_rms.write("ii=%d,jj=%d:\trms=%f\n"%(ii,jj,rms(new_data)))
        rms_mean += rms(new_data)
        count += 1

        if 0:
            plot(new_data)
            show()

hdulist[0].data = new_data
#synchronize header
hdulist[0].header.add_history(time.asctime(time.localtime())+'by python baseline_sin.py of GY')

hdulist.writeto(fname_o)

f_rms.write("mean rms:\t%f\n"%(rms_mean/count))
f_rms.close()

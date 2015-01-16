#!/usr/bin/env python
#filename=column density distribution
#2014/05/21
#---------------------------------------by Gao John

from __future__ import division
from pylab import *
import pyfits
import pywcsgrid2
import mpl_toolkits.axisartist as axisartist
from mpl_toolkits.axes_grid import make_axes_locatable
import matplotlib.colors as mcolors
import pywcs
from kapteyn import wcs

def locate_line(nx,ny,k,b,d_limit):
    if (b>=0 and b<ny-1):
        index_x = [0]
        index_y = [int(round(b))]
    elif (b<0):
        index_x = [int(round(-b/k))]
        index_y = [0]

    else:
        index_x = [int(round((ny-1-b)/k))]
        index_y = [ny-1]

    y_now = index_y[0]

    for i2 in xrange(nx):
        i = i2+index_x[0]  #x axis not begin at 0, but index_x[0]

        if (k>0 and y_now<ny-1):
            upper_run = True
            if (i2!=0): y_now-=1
        else:
            upper_run = False

        if (k<0 and y_now>0):
            down_run = True
            if (i2!=0): y_now+=1  #last y_now is not what I need
        else:
            down_run = False

        d = (1.+1./(k*k))*(k*i+b-y_now)*(k*i+b-y_now)  #distance^2
#        print k,b,d,d_limit
#        print i,y_now,upper_run,down_run
        if (i2==0 and d>2*d_limit):
            print 'warning initial d>d_limit'

#        if (d<2*d_limit):
#            index_x.append(i)
#            index_y.append(y_now)
#
#only one direction considering k
        if ((i!=0 and not (upper_run or down_run)) or (i>=nx)):
            break

#        print index_x,index_y
#        print upper_run,down_run
        running = True
        while running:
#            print i,y_now,k*i+b
            if (upper_run and y_now<ny-1):
                upper_run = True
            else:
                upper_run = False

            if (down_run and y_now>0):
                down_run = True
            else:
                down_run = False

            if (upper_run): 
                y_now+=1 
                d = (1.+1./(k*k))*(k*i+b-y_now)*(k*i+b-y_now)  #distance^2
                if (d>2*d_limit):
                    if (y_now>k*i+b):
                        upper_run = False
                    else:
                        pass

                else:
                    index_x.append(i)
                    index_y.append(y_now)


            if (down_run):
                y_now-=1
                d = (1.+1./(k*k))*(k*i+b-y_now)*(k*i+b-y_now)  #distance^2
                if (d>2*d_limit):
                    if (y_now<k*i+b):
                        down_run = False
                    else:
                        down_run = True

                else:
                    index_x.append(i)
                    index_y.append(y_now)


#            print upper_run, down_run
            running = upper_run or down_run

#    print array(index_y)
    return array(index_x),array(index_y)

#co_file===================================================
co_file = pyfits.open("kes41_12_-55--45.fits")
co_header = co_file[0].header
co_data = co_file[0].data
(ny, nx) = co_data.shape
co_wcs = pywcs.WCS(co_header)
print nx,ny

#pick x, y=========================================
proj = wcs.Projection(co_header)
pick_x_gal = 337.8      #galactic
pick_y_gal = -0.09
delt_sky_x = co_header['CDELT1']
delt_sky_y = co_header['CDELT2']
if (abs(co_header['CDELT1']) != abs(co_header['CDELT2'])):
    print 'warning: different between x and y'

trans = wcs.Transformation(wcs.galactic, proj.skysys)
(pick_x_co, pick_y_co) = trans.transform((pick_x_gal,pick_y_gal))
(pick_x, pick_y) = proj.topixel((pick_x_co, pick_y_co))
print pick_x,pick_y

(pick_x_co_tmp, pick_y_co_tmp) = trans.transform((pick_x_gal+5.*delt_sky_x,pick_y_gal))
(pick_x_tmp, pick_y_tmp) = proj.topixel((pick_x_co_tmp, pick_y_co_tmp))

if (pick_x_tmp == pick_x):
    print 'Error: vertical line on the map, you can get it a simpler way!'
    raise

k1 = (pick_y_tmp-pick_y)/(pick_x_tmp-pick_x)
b1 = pick_y-k1*pick_x

(pick_x_co_tmp, pick_y_co_tmp) = trans.transform((pick_x_gal,pick_y_gal+5.*delt_sky_y))
(pick_x_tmp, pick_y_tmp) = proj.topixel((pick_x_co_tmp, pick_y_co_tmp))
if (pick_x_tmp == pick_x):
    print 'Error: vertical line on the map, you can get it a simpler way!'
    raise

k2 = (pick_y_tmp-pick_y)/(pick_x_tmp-pick_x)
b2 = pick_y-k2*pick_x

#^lines of 2 different slope

d_limit = 1./4. #pixel

output1_x, output1_y = locate_line(nx,ny,k1,b1,d_limit)
output2_x, output2_y = locate_line(nx,ny,k2,b2,d_limit)

figure(0)
plot(output1_x,output1_y,'b-',output2_x,output2_y,'r-')
show()

output1 = zip(output1_x,output1_y)
output2 = zip(output2_x,output2_y)

output1_fk5 = proj.toworld(output1)
output2_fk5 = proj.toworld(output2)

trans_r = wcs.Transformation(proj.skysys,wcs.galactic)

output1_gal = trans_r.transform(output1_fk5)
output2_gal = trans_r.transform(output2_fk5)

output1_gal_ar = array(output1_gal)
output2_gal_ar = array(output2_gal)

print len(output1),len(output2)

output1_val = zeros(len(output1))
output2_val = zeros(len(output2))

for ind in xrange(len(output1)):
    output1_val[ind] = co_data.transpose()[output1[ind]]

for ind in xrange(len(output2)):
    output2_val[ind] = co_data.transpose()[output2[ind]]
#output2_val = co_data.transpose()[output2]
#plot(output1_gal_ar[:,0],output1_val,drawstyle='steps')

fig = figure(1)
fig.set_figheight(6.125)
fig.set_figwidth(16.25)
ax1 = fig.add_subplot(121)
ax1.plot(output1_gal_ar[:,0],output1_val*1.8e20/1e21,'b-',drawstyle='steps')
ax1.set_xlim(min(output1_gal_ar[:,0]),max(output1_gal_ar[:,0]))
#ax1.set_ylim(22*1.8e20,45*1.8e20)
ax1.set_xlabel('Galactic Latitude (deg)')
ax1.set_ylabel(r'N(H$_2$) (10$^{21}$ cm$^{-2}$)')
sp1 = ax1.axvspan(337.742247955,337.875792045,facecolor='c',alpha=0.3)
llist = [337.742,337.80,337.876]
llist_label = ['337.742','337.80','337.876']
ax1.set_xticks(llist)
ax1.set_xticklabels(llist_label)
ax2 = fig.add_subplot(122)
ax2.plot(output2_gal_ar[:,1],output2_val*1.8e20/1e21,'b-',drawstyle='steps')
ax2.set_xlim(min(output2_gal_ar[:,1]),max(output2_gal_ar[:,1]))
ax2.set_xlabel('Galactic Longitude (deg)')
ax2.set_ylabel(r'N(H$_2$) (10$^{21}$ cm$^{-2}$)')
#ax2.set_ylim(22*1.8e20,45*1.8e20)
sp2 = ax2.axvspan(-0.1405885325,-0.0492498775,facecolor='c',alpha=0.1)
blist = [-0.1406,-0.10,-0.0492]
blist_labels = [-0.1406,-0.10,-0.0492]
ax2.set_xticks(blist)
ax2.set_xticklabels(blist_labels)
fig.savefig('columndensity.eps')

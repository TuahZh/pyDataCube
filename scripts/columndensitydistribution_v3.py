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
co_file = pyfits.open("kes41_12_-70--40.fits")
co_header = co_file[0].header
co_data = co_file[0].data
(ny, nx) = co_data.shape
co_wcs = pywcs.WCS(co_header)
print nx,ny

#pick x, y=========================================
proj = wcs.Projection(co_header)
pick1_x1_fk5 = 249.825264
pick1_y1_fk5 = -46.907573
pick1_x2_fk5 = 249.6933
pick1_y2_fk5 = -47.00511
pick2_x1_fk5 = 249.72957
pick2_y1_fk5 = -46.928096
pick2_x2_fk5 = 249.815
pick2_y2_fk5 = -46.992042

delt_sky_x = co_header['CDELT1']
delt_sky_y = co_header['CDELT2']
if (abs(co_header['CDELT1']) != abs(co_header['CDELT2'])):
    print 'warning: different between x and y'

#trans = wcs.Transformation(wcs.galactic, proj.skysys)
#(pick_x_co, pick_y_co) = trans.transform((pick_x_gal,pick_y_gal))
(pick1_x1, pick1_y1) = proj.topixel((pick1_x1_fk5, pick1_y1_fk5))
(pick1_x2, pick1_y2) = proj.topixel((pick1_x2_fk5, pick1_y2_fk5))
(pick2_x1, pick2_y1) = proj.topixel((pick2_x1_fk5, pick2_y1_fk5))
(pick2_x2, pick2_y2) = proj.topixel((pick2_x2_fk5, pick2_y2_fk5))

#(pick_x_co_tmp, pick_y_co_tmp) = trans.transform((pick_x_gal+5.*delt_sky_x,pick_y_gal))
#(pick_x_tmp, pick_y_tmp) = proj.topixel((pick_x_co_tmp, pick_y_co_tmp))

#if (pick_x_tmp == pick_x):
#    print 'Error: vertical line on the map, you can get it a simpler way!'
#    raise

k1 = (pick1_y2-pick1_y1)/(pick1_x2-pick1_x1)
b1 = pick1_y2-k1*pick1_x2
k2 = (pick2_y2-pick2_y1)/(pick2_x2-pick2_x1)
b2 = pick2_y2-k2*pick2_x2

#(pick_x_co_tmp, pick_y_co_tmp) = trans.transform((pick_x_gal,pick_y_gal+5.*delt_sky_y))
#(pick_x_tmp, pick_y_tmp) = proj.topixel((pick_x_co_tmp, pick_y_co_tmp))
#if (pick_x_tmp == pick_x):
#    print 'Error: vertical line on the map, you can get it a simpler way!'
#    raise
#
#k2 = (pick_y_tmp-pick_y)/(pick_x_tmp-pick_x)
#b2 = pick_y-k2*pick_x
#
#^lines of 2 different slope

d_limit = 1./4. #pixel

output1_x, output1_y = locate_line(nx,ny,k1,b1,d_limit)
output2_x, output2_y = locate_line(nx,ny,k2,b2,d_limit)

figure(0)
plot(output1_x,output1_y,'b-',output2_x,output2_y,'r--')
show()

output1 = zip(output1_x,output1_y)
output2 = zip(output2_x,output2_y)

output1_fk5 = proj.toworld(output1)
output2_fk5 = proj.toworld(output2)
#trans_r = wcs.Transformation(proj.skysys,wcs.galactic)
#
#output1_gal = trans_r.transform(output1_fk5)
#output2_gal = trans_r.transform(output2_fk5)
#
#output1_gal_ar = array(output1_gal)
#output2_gal_ar = array(output2_gal)
#

print len(output1)

output1_val = zeros(len(output1))
output2_val = zeros(len(output2))

for ind in xrange(len(output1)):
    output1_val[ind] = co_data.transpose()[output1[ind]]

for ind in xrange(len(output2)):
    output2_val[ind] = co_data.transpose()[output2[ind]]
#output2_val = co_data.transpose()[output2]
#plot(output1_gal_ar[:,0],output1_val,drawstyle='steps')

#angle from center
if k1==k2:
    print "k1==k2!"
    raise

center_of_line = proj.toworld(((b2-b1)/(k1-k2),(k1*b2-k2*b1)/(k1-k2)))
print center_of_line
output1_angle = zeros(len(output1))
output2_angle = zeros(len(output2))

for ind in xrange(len(output1)):
#    output_angle[ind] = sqrt((output_fk5[ind][0]-center_of_line[0])* \
#            (output_fk5[ind][0]-center_of_line[0])+ \
#            (output_fk5[ind][1]-center_of_line[1])*(output_fk5[ind][1]-center_of_line[1]))* \
#            (output_fk5[ind][0]-center_of_line[0])
    output1_angle[ind] = linalg.norm(array(output1_fk5[ind])-\
#distance:linalg.norm
            array(center_of_line))*\
            (-sign(output1_fk5[ind][0]-center_of_line[0]))

for ind in xrange(len(output2)):
    output2_angle[ind] = linalg.norm(array(output2_fk5[ind])-\
            array(center_of_line))*\
            (-sign(output2_fk5[ind][0]-center_of_line[0]))

fig = figure(1)
fig.set_figheight(6.125)
fig.set_figwidth(16.25)
ax1 = fig.add_subplot(121,adjustable='box')
ax1.plot(output1_angle*60.,output1_val/0.42*1.8e20/1e21,'b-',drawstyle='steps')
ax1.set_xlim(min(output1_angle*60.),max(output1_angle*60.))
#ax1.set_ylim(22*1.8e20,45*1.8e20)
ax1.set_xlabel('Offset from the reference point (arcmin)',fontsize='x-large')
ax1.set_ylabel(r'N(H$_2$) (10$^{21}$ cm$^{-2}$)',fontsize='x-large')
#angle1 = 60.*sqrt((pick1_x1_fk5-center_of_line[0])*(pick1_x1_fk5-center_of_line[0])+\
#        (pick1_y1_fk5-center_of_line[1])*(pick1_y1_fk5-center_of_line[1]))
angle1 = 60.*linalg.norm(array([pick1_x1_fk5,pick1_y1_fk5])-array(center_of_line))*sign(pick1_x1_fk5-center_of_line[0])
#in unit of minute
angle2 = 60.*linalg.norm(array([pick1_x2_fk5,pick1_y2_fk5])-array(center_of_line))*sign(pick1_x2_fk5-center_of_line[0])
#angle2 = -60.*sqrt((pick_x2_fk5-center_of_line[0])*(pick_x2_fk5-center_of_line[0])+\
#        (pick_y2_fk5-center_of_line[1])*(pick_y2_fk5-center_of_line[1]))
sp1 = ax1.axvspan(angle1,angle2,facecolor='#d5d5d5',edgecolor='w',alpha=0.3)
xlist = [-6.0,-3.0,0.0,3.0,6.0]
xlist_label = ['-6.0','-3.0','0.0','3.0','6.0']
ax1.set_xticks(xlist)
ax1.set_xticklabels(xlist_label,size='large')
ylist = [28.0,32.0,36.0,40.0,44.0]
ylist_label = ['28.0','32.0','36.0','40.0','44.0']
ax1.set_yticks(ylist)
ax1.set_yticklabels(ylist_label,size='large')
text1="NE-SW"
ax1.text(-6,41.5,text1,fontsize=14,ha='left')

ax2 = fig.add_subplot(122,adjustable='box')
ax2.plot(output2_angle*60.,output2_val/0.42*1.8e20/1e21,'b-',drawstyle='steps')
ax2.set_xlim(min(output2_angle*60.),max(output2_angle*60.))
#ax1.set_ylim(22*1.8e20,45*1.8e20)
ax2.set_xlabel('Offset from the reference point (arcmin)',fontsize='x-large')
#ax2.set_ylabel(r'N(H$_2$) (10$^{21}$ cm$^{-2}$)',fontsize='x-large')
#angle1 = 60.*sqrt((pick1_x1_fk5-center_of_line[0])*(pick1_x1_fk5-center_of_line[0])+\
#        (pick1_y1_fk5-center_of_line[1])*(pick1_y1_fk5-center_of_line[1]))
angle1 = 60.*linalg.norm(array([pick2_x1_fk5,pick2_y1_fk5])-array(center_of_line))*sign(pick2_x1_fk5-center_of_line[0])
#in unit of minute
angle2 = 60.*linalg.norm(array([pick2_x2_fk5,pick2_y2_fk5])-array(center_of_line))*sign(pick2_x2_fk5-center_of_line[0])
#angle2 = -60.*sqrt((pick_x2_fk5-center_of_line[0])*(pick_x2_fk5-center_of_line[0])+\
#        (pick_y2_fk5-center_of_line[1])*(pick_y2_fk5-center_of_line[1]))
sp2 = ax2.axvspan(angle1,angle2,facecolor='#d5d5d5',edgecolor='w',alpha=0.3)
xlist = [-6.0,-3.0,0.0,3.0,6.0]
xlist_label = ['-6.0','-3.0','0.0','3.0','6.0']
ax2.set_xticks(xlist)
ax2.set_xticklabels(xlist_label,size='large')
ylist = [15.0,20.0,25.0,30.0,35.0,40.0,45.0]
ylist_label = ['15.0','20.0','25.0','30.0','35.0','40.0','45.0']
ax2.set_yticks(ylist)
ax2.set_yticklabels(ylist_label,size='large')
text2='NW-SE'
ax2.text(-6.5,40,text2,fontsize=14,ha='left')

fig.savefig('columndensity_v3.eps')
#show()

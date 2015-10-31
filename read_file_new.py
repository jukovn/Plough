# -*- coding: utf-8 -*-
from numpy import array;
from numpy import zeros;
from numpy import pi;
from numpy import sin;
from numpy import cos;
from numpy import sqrt;
from numpy import copy;
from numpy import power;
from numpy import e;
import matplotlib.pyplot as plt


num_w=210;
#num_w=1991;
num_extr=100;

w_ar1=zeros(num_w);
ampl_ar1=zeros((num_w*num_extr));
i_w=0;
i=0;
#file_out=open('C:\\Domashnie_zadania\\Freza_ploskoe\\Gotovo\\dyn_koleb_0.3-2,3,0.01.txt','r');
file_out=open('C:\\Domashnie_zadania\\Freza_ploskoe\\Program\\Results_1\\without_fl_x.txt','r');

for line in file_out:
    if i_w<num_w:
        if (i%(num_extr+1))==0:
            w_ar1[i_w]=float(line);
            i_w=i_w+1;
            print i_w,num_w-10,line;
        else:
            ampl_ar1[(i_w-1)*num_extr+i%(num_extr+1)-1]=float(line);
        i=i+1;
file_out.close();

w_pl1=zeros((num_w*num_extr));
for i in range(num_w*num_extr):
    w_pl1[i]=w_ar1[i/100];


for i in range(num_w*num_extr):
    if ampl_ar1[i]>=1:
        ampl_ar1[i]=1;
    if ampl_ar1[i]<=-1:
        ampl_ar1[i]=-1;
        
        
###
###
###




w_ar2=zeros(num_w);
ampl_ar2=zeros((num_w*num_extr));
i_w=0;
i=0;
file_out=open('C:\\Domashnie_zadania\\Freza_ploskoe\\Program\\Results_1\\eps_h=0.5_x.txt','r');
for line in file_out:
    if i_w<num_w:
        if (i%(num_extr+1))==0:
            w_ar2[i_w]=float(line);
            i_w=i_w+1;
        else:
            ampl_ar2[(i_w-1)*num_extr+i%(num_extr+1)-1]=float(line);
        i=i+1;
file_out.close();

w_pl2=zeros((num_w*num_extr));
for i in range(num_w*num_extr):
    w_pl2[i]=w_ar2[i/100];

for i in range(num_w*num_extr):
    if ampl_ar2[i]>=1:
        ampl_ar2[i]=1;
    if ampl_ar2[i]<=-1:
        ampl_ar2[i]=-1;



###
###
###



w_ar3=zeros(num_w);
ampl_ar3=zeros((num_w*num_extr));
i_w=0;
i=0;
file_out=open('C:\\Domashnie_zadania\\Freza_ploskoe\\Program\\Results_1\\eps_h=0.2_x.txt','r');
for line in file_out:
    if i_w<num_w:
        if (i%(num_extr+1))==0:
            w_ar3[i_w]=float(line);
            i_w=i_w+1;
        else:
            ampl_ar3[(i_w-1)*num_extr+i%(num_extr+1)-1]=float(line);
        i=i+1;
file_out.close();

w_pl3=zeros((num_w*num_extr));
for i in range(num_w*num_extr):
    w_pl3[i]=w_ar3[i/100];

for i in range(num_w*num_extr):
    if ampl_ar3[i]>=1:
        ampl_ar3[i]=1;
    if ampl_ar3[i]<=-1:
        ampl_ar3[i]=-1;


"""

w_ar4=zeros(num_w);
ampl_ar4=zeros((num_w*num_extr));
i_w=0;
i=0;
file_out=open('C:\\Domashnie_zadania\\Freza_ploskoe\\zadnia_gran\\fr_out_y-2.txt','r');
for line in file_out:
    if (i%(num_extr+1))==0:
        w_ar4[i_w]=float(line);
        i_w=i_w+1;
    else:
        ampl_ar4[(i_w-1)*num_extr+i%(num_extr+1)-1]=float(line);
    i=i+1;
file_out.close();

w_pl4=zeros((num_w*num_extr));
for i in range(num_w*num_extr):
    w_pl4[i]=w_ar4[i/100];


for i in range(num_w*num_extr):
    if ampl_ar4[i]>=1:
        ampl_ar4[i]=1;
    if ampl_ar4[i]<=-1:
        ampl_ar4[i]=-1;












"""


###
###
###


plt.figure(2000);
plt.subplot(1,2,1);
plt.plot(w_pl3,ampl_ar3, 'ro',markerfacecolor='r');
plt.xlabel('w_otn')
plt.ylabel('x_dyn')
plt.title('fl_otn=0.2')

plt.subplot(1,2,2);
plt.plot(w_pl1,ampl_ar1, 'ro',markerfacecolor='r');
plt.xlabel('w_otn')
plt.ylabel('x_dyn')
plt.title('fl_otn=0.5')
"""
plt.subplot(1,3,3);
plt.plot(w_pl3,ampl_ar3, 'ro',markerfacecolor='r');
plt.xlabel('w_otn')
plt.ylabel('x_dyn')
plt.title('fl_otn=0.2')

"""
"""
plt.figure(2001);
plt.subplot(1,2,2);
plt.plot(w_pl2,ampl_ar2, 'ro',markerfacecolor='r');
plt.xlabel('w_otn')
plt.ylabel('h_fl')
plt.title('Nadpis')

plt.figure(2002);
plt.subplot(1,2,2);
plt.plot(w_pl3,ampl_ar3, 'ro',markerfacecolor='r');
plt.xlabel('w_otn')
plt.ylabel('h')
plt.title('Nadpis')

plt.figure(2003);
plt.subplot(1,2,2);
plt.plot(w_pl4,ampl_ar4, 'ro',markerfacecolor='r');
plt.xlabel('w_otn')
plt.ylabel('y_dyn')
plt.title('Nadpis')
"""
plt.show();
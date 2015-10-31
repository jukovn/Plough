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
            #print i_w,num_w-10,line;
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

new_ampl_ar1=zeros((num_w*num_extr));
for i in range(num_w):
    new_ampl_ar1[(0+100*i):(50+100*i)]=max(ampl_ar1[(100*i):(100*(i+1))]);
    new_ampl_ar1[(50+100*i):(100+100*i)]=min(ampl_ar1[(100*i):(100*(i+1))]);
    
    
        
###
###
###




w_ar2=zeros(num_w);
ampl_ar2=zeros((num_w*num_extr));
i_w=0;
i=0;
file_out=open('C:\\Domashnie_zadania\\Freza_ploskoe\\Program\\Results_1\\eps_h=0.8_x.txt','r');
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

new_ampl_ar2=zeros((num_w*num_extr));
for i in range(num_w):
    new_ampl_ar2[(0+100*i):(50+100*i)]=max(ampl_ar2[(100*i):(100*(i+1))]);
    new_ampl_ar2[(50+100*i):(100+100*i)]=min(ampl_ar2[(100*i):(100*(i+1))]);

###
###
###



w_ar3=zeros(num_w);
ampl_ar3=zeros((num_w*num_extr));
i_w=0;
i=0;
file_out=open('C:\\Domashnie_zadania\\Freza_ploskoe\\Program\\Results_1\\eps_h=0.6_x.txt','r');
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


new_ampl_ar3=zeros((num_w*num_extr));
for i in range(num_w):
    new_ampl_ar3[(0+100*i):(50+100*i)]=max(ampl_ar3[(100*i):(100*(i+1))]);
    new_ampl_ar3[(50+100*i):(100+100*i)]=min(ampl_ar3[(100*i):(100*(i+1))]);




w_ar4=zeros(num_w);
ampl_ar4=zeros((num_w*num_extr));
i_w=0;
i=0;
file_out=open('C:\\Domashnie_zadania\\Freza_ploskoe\\Program\\Results_1\\eps_h=0.4_x.txt','r');
for line in file_out:
    if i_w<num_w:
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

new_ampl_ar4=zeros((num_w*num_extr));
for i in range(num_w):
    new_ampl_ar4[(0+100*i):(50+100*i)]=max(ampl_ar4[(100*i):(100*(i+1))]);
    new_ampl_ar4[(50+100*i):(100+100*i)]=min(ampl_ar4[(100*i):(100*(i+1))]);






w_ar5=zeros(num_w);
ampl_ar5=zeros((num_w*num_extr));
i_w=0;
i=0;
file_out=open('C:\\Domashnie_zadania\\Freza_ploskoe\\Program\\Results_1\\eps_h=0.2_x.txt','r');
for line in file_out:
    if i_w<num_w:
        if (i%(num_extr+1))==0:
            w_ar5[i_w]=float(line);
            i_w=i_w+1;
        else:
            ampl_ar5[(i_w-1)*num_extr+i%(num_extr+1)-1]=float(line);
        i=i+1;
file_out.close();

w_pl5=zeros((num_w*num_extr));
for i in range(num_w*num_extr):
    w_pl5[i]=w_ar5[i/100];


for i in range(num_w*num_extr):
    if ampl_ar5[i]>=1:
        ampl_ar5[i]=1;
    if ampl_ar5[i]<=-1:
        ampl_ar5[i]=-1;

new_ampl_ar5=zeros((num_w*num_extr));
for i in range(num_w):
    new_ampl_ar5[(0+100*i):(50+100*i)]=max(ampl_ar5[(100*i):(100*(i+1))]);
    new_ampl_ar5[(50+100*i):(100+100*i)]=min(ampl_ar5[(100*i):(100*(i+1))]);










###
###
###


plt.figure(2000);

plt.plot(w_pl1,new_ampl_ar1, 'ro',markerfacecolor='b');
plt.plot(w_pl2,new_ampl_ar2, 'ro',markerfacecolor='g');
plt.plot(w_pl3,new_ampl_ar3, 'ro',markerfacecolor='r');
plt.plot(w_pl4,new_ampl_ar4, 'ro',markerfacecolor='c');
plt.plot(w_pl5,new_ampl_ar5, 'ro',markerfacecolor='m');





plt.xlabel('w_otn')
plt.ylabel('x_dyn')
plt.title('fl_otn=0.5')







"""
plt.subplot(1,2,2);
plt.plot(w_pl1,new_ampl_ar1, 'ro',markerfacecolor='r');
plt.xlabel('w_otn')
plt.ylabel('x_dyn')
plt.title('fl_otn=0.2')
"""
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
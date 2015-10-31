# -*- coding: utf-8 -*-

def k_xy(x1,y1,x2,y2):
    #if (x1==x2):
    #    print "och ploha 1",x1,y1,x2,y2;
    return (y2-y1)/(x2-x1);

def b_xy(x1,y1,x2,y2):
    return (y1*x2-x1*y2)/(x2-x1);

def x_per(k1,b1,k2,b2):
    #if (k2==k1):
    #    print "och ploha 2";
    return (b1-b2)/(k2-k1);

def y_x(k,b,x):
    return (k*x+b);

def h_rad(rez_x,rez_y,x_p,y_p):
    return sqrt((rez_x-x_p)**2+(rez_y-y_p)**2);

def h_det(z,dx,rez_x,rez_y,xc,yc):
    zi=int(rez_x/dx);
    k_s=k_xy(zi*dx,z[zi],(zi+1)*dx,z[(zi+1)]);
    b_s=b_xy(zi*dx,z[zi],(zi+1)*dx,z[(zi+1)]);
    z_s=y_x(k_s,b_s,rez_x);
            
    if (z_s<=rez_y):
        #если текущая точка лежит вне поверхности
        h=0;
                                                            
    else:
         #если текушая точка лежит внутри поверхности
            
        k_f=k_xy(xc,yc,rez_x,rez_y);
        b_f=b_xy(xc,yc,rez_x,rez_y);
        count_p=0;
        j=0;
        if (rez_x<=xc):
            while (count_p==0):
                #print "1",zi+j,rez_x,xc;
                k_s=k_xy((zi+j)*dx,z[(zi+j)],(zi+j+1)*dx,z[(zi+j+1)]);
                b_s=b_xy((zi+j)*dx,z[(zi+j)],(zi+j+1)*dx,z[(zi+j+1)]);
                x_p=x_per(k_f,b_f,k_s,b_s);
                if ((x_p>=(zi+j)*dx) and (x_p<=(zi+j+1)*dx)):
                    count_p=1;
                else:
                    j=j+1;
        else:
            while (count_p==0):
                #print "2",zi+1-j;
                k_s=k_xy((zi+1-j)*dx,z[(zi+1-j)],(zi-j)*dx,z[(zi-j)]);
                b_s=b_xy((zi+1-j)*dx,z[(zi+1-j)],(zi-j)*dx,z[(zi-j)]);
                x_p=x_per(k_f,b_f,k_s,b_s);
                if ((x_p>=(zi-j)*dx) and (x_p<=(zi-j+1)*dx)):
                    count_p=1;
                else:
                    j=j+1;
                
        h=h_rad(rez_x,rez_y,x_p,y_x(k_f,b_f,x_p));
    return h;

def dyn_per(a,b,x_0,V_0,n,p1,m,dt):
    #print b,a*dt+b;
    #print a,b,x_0,V_0,n,p1,m,dt;
    #F=a+b*t
    #x_0 поступает в мм
    #V_0 поступает в мм\с
    #перевод в м и м\с
    #x_0=x_0/1000;
    #V_0=V_0/1000;
    out_ar=zeros((2));
    
    #o_o=power(e,-n*dt)*(x_0*cos(p1*dt)+(V_0+n*x_0)/p1*sin(p1*dt));
    #c_n=b/(m*p1*(n**2+p1**2))*(p1-power(e,-n*dt)*(p1*cos(p1*dt)+n*sin(p1*dt)));
    #out_ar[0]=o_o+c_n;
    
    o_o=power(e,-n*dt)*(x_0*cos(p1*dt)+(V_0+n*x_0)/p1*sin(p1*dt));
    c_n_1=-b/(n**2+p1**2)*(-power(e,n*dt)*p1+p1*cos(p1*dt)+n*sin(p1*dt));
    c_n_2=a/power(n**2+p1**2,2)*(power(e,n*dt)*p1*(-2*n+(n**2+p1**2)*dt)+2*n*p1*cos(p1*dt)+(n**2-p1**2)*sin(p1*dt));
    out_ar[0]=1.0*(o_o+power(e,-n*dt)/(m*p1)*(c_n_1+c_n_2));
    
    #v_c_n=b*sin(p1*dt)*power(e,-n*dt)/(m*p1);
    #v_o_o=power(e,-n*dt)*(V_0*(cos(p1*dt)-n/p1*sin(p1*dt))-x_0*(n**2/p1*sin(p1*dt)+p1*sin(p1*dt)));
    #out_ar[1]=v_o_o+v_c_n;
    
    
    cnV_1=p1*1/power(n**2+p1**2,2)*(power(e,n*dt)*(b*n*(n**2+p1**2)+a*(-n**2+p1**2+n*(n**2+p1**2)*dt))+(n**2*(a-b*n)-(a+b*n)*p1**2)*cos(p1*dt)+p1*(-2*a*n+b*(n**2+p1**2))*sin(p1*dt));
    V_k=-n*(o_o+power(e,-n*dt)/(m*p1)*(c_n_1+c_n_2))+power(e,-n*dt)*(-x_0*p1*sin(p1*dt)+(V_0+n*x_0)*cos(p1*dt))+power(e,-n*dt)/(m*p1)*cnV_1;
    out_ar[1]=1.0*V_k;
    
    return out_ar;
 
def a_F(F1,F2,dt):
    return (F2-F1)/dt;
    
def x_norm(k1,b1,x0,y0):
    return ((x0+k1*y0-k1*b1)/(k1**2+1));

def y_norm(k1,b1,x0,y0):
    return ((k1*x0+k1**2*y0+b1)/(k1**2+1));
    
def h_fl_det(rez_fin,rez_flank1,data_s1,xc,yc):
    
    data_s=data_s1-1;
    rez_flank=zeros((3,data_s));
    rez_flank[0]=rez_flank1[0][1:data_s1];
    rez_flank[1]=rez_flank1[1][1:data_s1];
    rez_flank[2]=rez_flank1[2][1:data_s1];
    
    
    h_fl_max=0;
    k_fl=k_xy(rez_fin[4],rez_fin[5],rez_fin[2],rez_fin[3]);
    b_fl=b_xy(rez_fin[4],rez_fin[5],rez_fin[2],rez_fin[3]);
    x_fl_1=min(rez_fin[4],rez_fin[2]);
    x_fl_2=max(rez_fin[4],rez_fin[2]);
    if (rez_fin[4]==x_fl_1):
        y_fl_1=rez_fin[5];
        y_fl_2=rez_fin[3];
    elif (rez_fin[4]==x_fl_2):
        y_fl_1=rez_fin[3];
        y_fl_2=rez_fin[5];
    
    
    
    mask=rez_flank[1]>y_x(k_fl,b_fl,rez_flank[0]);#проверка условия (2)
    x_osn_h=x_norm(k_fl,b_fl,rez_flank[0][mask],rez_flank[1][mask]);
    y_osn_h=y_norm(k_fl,b_fl,rez_flank[0][mask],rez_flank[1][mask]);
    mask1=x_osn_h>=x_fl_1;#проверка условия (3)
    mask2=x_osn_h[mask1]<=x_fl_2;#проверка условия (3)
    h_norm=sqrt((rez_flank[0][mask][mask1][mask2]-x_osn_h[mask1][mask2])**2+(rez_flank[1][mask][mask1][mask2]-y_osn_h[mask1][mask2])**2);
    if h_norm.size!=0:
        h_fl_max=max(h_norm);#проверка условия (4)
    else:
        h_fl_max=0;
    
    
    #определение точки пересечения линии передней грани и поверхности
    
    k_f=k_xy(rez_fin[2],rez_fin[3],xc,yc);
    b_f=b_xy(rez_fin[2],rez_fin[3],xc,yc);
    count_c=0;
    mask=rez_flank[2][0:data_s-1]//(2*pi)==rez_flank[2][1:data_s]//(2*pi);#проверка того, что алгоритм не "соединит" точки с разными количествами оборотов
    k_s_fl=k_xy(rez_flank[0][0:data_s-1][mask],rez_flank[1][0:data_s-1][mask],rez_flank[0][1:data_s][mask],rez_flank[1][1:data_s][mask]);
    b_s_fl=b_xy(rez_flank[0][0:data_s-1][mask],rez_flank[1][0:data_s-1][mask],rez_flank[0][1:data_s][mask],rez_flank[1][1:data_s][mask]);
    mask1=rez_fin[3]<=y_x(k_s_fl,b_s_fl,rez_fin[2]);#проверка, что вершина резца лежит "внутри" поверхности
    x_p=x_per(k_f,b_f,k_s_fl[mask1],b_s_fl[mask1]);
    mask2=x_p>=minimum(rez_flank[0][0:data_s-1][mask][mask1],rez_flank[0][1:data_s][mask][mask1]);
    mask3=x_p[mask2]<=maximum(rez_flank[0][0:data_s-1][mask][mask1][mask2],rez_flank[0][1:data_s][mask][mask1][mask2]);
    
    ###
    mask4=x_p[mask2][mask3]>=min(xc,rez_fin[2]);
    mask5=x_p[mask2][mask3][mask4]<=max(xc,rez_fin[2]);
    
    ###
    
    
    if x_p[mask2][mask3][mask4][mask5].size==1:
        k_s_per=k_s_fl[mask1][mask2][mask3][mask4][mask5][0];
        b_s_per=b_s_fl[mask1][mask2][mask3][mask4][mask5][0];
        x_per_fl=x_p[mask2][mask3][mask4][mask5][0];
        y_per_fl=y_x(k_s_per,b_s_per,x_per_fl);
        x_p_osn=x_norm(k_fl,b_fl,x_per_fl,y_per_fl);
        y_p_osn=y_norm(k_fl,b_fl,x_per_fl,y_per_fl);
        h_per=sqrt((x_p_osn-x_per_fl)**2+(y_p_osn-y_per_fl)**2);
        h_fl_max=max(h_fl_max,h_per);
    

                
      
    
    return h_fl_max;



def K_fl_cont(S,Krc):
    if S==0:
        return 0;
    else:
        #return 0.01*Krc/(1.8*10**(-7)+9.62*10**(-2)*S+2.11*10**2*S**2);
        #return 0.01*Krc/(1.8*10**(-7)+9.62*10**(-2)*S);
        #return 0.1*Krc/(10**-5+10*S);
        #return 0.1*Krc/(10**-5+10**4*S**2);
        return 0.1*Krc/(S);
        
def K_fl_cont1(S,S0,Krc):
    if S==0:
        return 0;
    else:
        return Krc/(S0-S);


def s_fl_det(rez_fin,rez_flank1,data_s1,xc,yc):
    
    data_s=data_s1-1;
    rez_flank=zeros((3,data_s));
    rez_flank[0]=rez_flank1[0][1:data_s1];
    rez_flank[1]=rez_flank1[1][1:data_s1];
    rez_flank[2]=rez_flank1[2][1:data_s1];
    
    x_per_r=0;
    x_per_l=0;
    
    s_fl=0;
    
    k_fl=k_xy(rez_fin[4],rez_fin[5],rez_fin[2],rez_fin[3]);
    b_fl=b_xy(rez_fin[4],rez_fin[5],rez_fin[2],rez_fin[3]);
    
    #mask=rez_flank[2][0:data_s-1]//(2*pi)==rez_flank[2][1:data_s]//(2*pi);#проверка того, что алгоритм не "соединит" точки с разными количествами оборотов
    mask=(rez_flank[2]//(2*pi))==(rez_flank[2][0]//(2*pi));#рассматриваем только точки с тем же количеством полных оборотов, что и у текущей точки
    #print mask;
    k_s_fl=k_xy(rez_flank[0][mask][:-1],rez_flank[1][mask][:-1],rez_flank[0][mask][1:],rez_flank[1][mask][1:]);
    b_s_fl=b_xy(rez_flank[0][mask][:-1],rez_flank[1][mask][:-1],rez_flank[0][mask][1:],rez_flank[1][mask][1:]);
    #k_s_fl=k_xy(rez_flank[0][0:data_s-1][mask],rez_flank[1][0:data_s-1][mask],rez_flank[0][1:data_s][mask],rez_flank[1][1:data_s][mask]);
    #b_s_fl=b_xy(rez_flank[0][0:data_s-1][mask],rez_flank[1][0:data_s-1][mask],rez_flank[0][1:data_s][mask],rez_flank[1][1:data_s][mask]);
    
    
    #определение точки пересечения линии передней грани и поверхности
    
    k_f=k_xy(rez_fin[2],rez_fin[3],xc,yc);
    b_f=b_xy(rez_fin[2],rez_fin[3],xc,yc);
    
    mask1=rez_fin[3]<=y_x(k_s_fl,b_s_fl,rez_fin[2]);#проверка, что вершина резца лежит "внутри" поверхности
    x_p=x_per(k_f,b_f,k_s_fl[mask1],b_s_fl[mask1]);
    mask2=x_p>=minimum(rez_flank[0][mask][:-1][mask1],rez_flank[0][mask][1:][mask1]);
    mask3=x_p[mask2]<=maximum(rez_flank[0][mask][:-1][mask1][mask2],rez_flank[0][mask][1:][mask1][mask2]);
    mask4=x_p[mask2][mask3]>=min(xc,rez_fin[2]);
    mask5=x_p[mask2][mask3][mask4]<=max(xc,rez_fin[2]);
    
 
    if x_p[mask2][mask3][mask4][mask5].size!=0:
        k_s_r=k_s_fl[mask1][mask2][mask3][mask4][mask5][0];
        b_s_r=b_s_fl[mask1][mask2][mask3][mask4][mask5][0];
        x_per_r=x_p[mask2][mask3][mask4][mask5][0];
        y_per_r=y_x(k_s_r,b_s_r,x_per_r);
        #print "1",x_per_r,y_per_r;
    else:
        x_p1=x_per(k_fl,b_fl,k_s_fl,b_s_fl);
        mask1=x_p1>=minimum(rez_flank[0][mask][:-1],rez_flank[0][mask][1:]);
        mask2=x_p1[mask1]<=maximum(rez_flank[0][mask][:-1][mask1],rez_flank[0][mask][1:][mask1]);
        
        mask3=x_p1[mask1][mask2]>=min(rez_fin[2],rez_fin[4]);
        mask4=x_p1[mask1][mask2][mask3]<=max(rez_fin[2],rez_fin[4]);
        
        if x_p1[mask1][mask2][mask3][mask4].size!=0:
            k_s_r=k_s_fl[mask1][mask2][mask3][mask4][0];
            b_s_r=b_s_fl[mask1][mask2][mask3][mask4][0];
            x_per_r=x_p1[mask1][mask2][mask3][mask4][0];
            y_per_r=y_x(k_s_r,b_s_r,x_per_r);
            #print "2",x_per_r,y_per_r;
    
    #определение точки пересечения "второй линии" режущей кромки и поверхности
    
    k_fb=k_xy(rez_fin[4],rez_fin[5],xc,yc);
    b_fb=b_xy(rez_fin[4],rez_fin[5],xc,yc);
    
    mask1=rez_fin[5]<=y_x(k_s_fl,b_s_fl,rez_fin[4]);#проверка, что "вторая" точка задней грани лежит "внутри" поверхности
    x_p=x_per(k_fb,b_fb,k_s_fl[mask1],b_s_fl[mask1]);
    mask2=x_p>=minimum(rez_flank[0][mask][:-1][mask1],rez_flank[0][mask][1:][mask1]);
    mask3=x_p[mask2]<=maximum(rez_flank[0][mask][:-1][mask1][mask2],rez_flank[0][mask][1:][mask1][mask2]);
    mask4=x_p[mask2][mask3]>=min(xc,rez_fin[4]);
    mask5=x_p[mask2][mask3][mask4]<=max(xc,rez_fin[4]);
    
    
    if x_p[mask2][mask3][mask4][mask5].size!=0:
        k_s_l=k_s_fl[mask1][mask2][mask3][mask4][mask5][0];
        b_s_l=b_s_fl[mask1][mask2][mask3][mask4][mask5][0];
        x_per_l=x_p[mask2][mask3][mask4][mask5][0];
        y_per_l=y_x(k_s_l,b_s_l,x_per_l);
        #print "3",x_per_l,y_per_l;
    else:
        x_p1=x_per(k_fl,b_fl,k_s_fl,b_s_fl);
        mask1=x_p1>=minimum(rez_flank[0][mask][:-1],rez_flank[0][mask][1:]);
        mask2=x_p1[mask1]<=maximum(rez_flank[0][mask][:-1][mask1],rez_flank[0][mask][1:][mask1]);
        
        mask3=x_p1[mask1][mask2]>=min(rez_fin[2],rez_fin[4]);
        mask4=x_p1[mask1][mask2][mask3]<=max(rez_fin[2],rez_fin[4]);
        
        if x_p1[mask1][mask2][mask3][mask4].size!=0:
            k_s_l=k_s_fl[mask1][mask2][mask3][mask4][-1];
            b_s_l=b_s_fl[mask1][mask2][mask3][mask4][-1];
            x_per_l=x_p1[mask1][mask2][mask3][mask4][-1];
            y_per_l=y_x(k_s_l,b_s_l,x_per_l);
            #print "4",x_per_l,y_per_l;
    
    if (x_per_r==0) or (x_per_l==0):
        s_fl=0;
    else:
        
        mask1=rez_flank[0][mask]<=x_per_r;
        mask2=rez_flank[0][mask][mask1]>=x_per_l;
        #print "popali",rez_flank[0][mask][mask1][mask2];
        
        s_fl_num=rez_flank[0][mask][mask1][mask2].size;
        #print "num",s_fl_num;
        
        if s_fl_num!=0:
            #подсчет "малой правой" площади
            s_fl=s_fl+0.5*(rez_fin[2]-x_per_r)*(y_per_r-y_x(k_fl,b_fl,x_per_r));#площадь малого правого треугольника
            x_s_r=rez_flank[0][mask][rez_flank[0][mask]<=x_per_r][0];#координаты бижайшей левой точки на поверности 
            y_s_r=rez_flank[1][mask][rez_flank[0][mask]<=x_per_r][0];#к вершине резца
            #print x_s_r,y_s_r;
            s_fl=s_fl+0.5*(x_per_r-x_s_r)*((y_per_r-y_x(k_fl,b_fl,x_per_r))+(y_s_r-y_x(k_fl,b_fl,x_s_r)));#площадь малой правой трапеции
            #print s_fl;
            
            #подсчет "малой левой" площади
            s_fl=s_fl+0.5*(x_per_l-rez_fin[4])*(y_per_l-y_x(k_fl,b_fl,x_per_l));#площадь малого левого треугольника
            x_s_l=rez_flank[0][mask][rez_flank[0][mask]>=x_per_l][-1];#координаты ближайшей правой точки на поверхности
            y_s_l=rez_flank[1][mask][rez_flank[0][mask]>=x_per_l][-1];#к "второй вершине" задней грани режущей кромки
            #print x_s_l,y_s_l;
            s_fl=s_fl+0.5*(x_s_l-x_per_l)*((y_per_l-y_x(k_fl,b_fl,x_per_l))+(y_s_l-y_x(k_fl,b_fl,x_s_l)));#площадь малой левой трапеции
            #print s_fl;
            
            s_fl_intern=0.5*(rez_flank[0][mask][mask1][mask2][:-1]-rez_flank[0][mask][mask1][mask2][1:])*((rez_flank[1][mask][mask1][mask2][1:]-y_x(k_fl,b_fl,rez_flank[0][mask][mask1][mask2][1:]))+(rez_flank[1][mask][mask1][mask2][:-1]-y_x(k_fl,b_fl,rez_flank[0][mask][mask1][mask2][:-1])));
            
            #print "6",s_fl_intern;    
            s_fl=s_fl+sum(s_fl_intern);
            
        else:
            s_fl=s_fl+0.5*(rez_fin[2]-x_per_r)*(y_per_r-y_x(k_fl,b_fl,x_per_r));#площадь правого треугольника
            
            s_fl=s_fl+0.5*(x_per_l-rez_fin[4])*(y_per_l-y_x(k_fl,b_fl,x_per_l));#площадь левого треугольника
            
            s_fl=s_fl+0.5*(x_per_r-x_per_l)*((y_per_r-y_x(k_fl,b_fl,x_per_r))+(y_per_l-y_x(k_fl,b_fl,x_per_l)));
            
            #print "10",s_fl;
    
    if s_fl<0:
        s_fl=0;
    
   
    return s_fl;

      
from numpy import array;
from numpy import zeros;
from numpy import pi;
from numpy import sin;
from numpy import cos;
from numpy import sqrt;
from numpy import copy;
from numpy import power;
from numpy import e;
from numpy import arcsin;
from numpy import insert;
from numpy import sort;
from numpy import argsort;
from numpy import minimum;
from numpy import maximum;
from numpy import linspace
import matplotlib.pyplot as plt;
import matplotlib.lines as mlines;


#задание начальных данных
L=150.0;#в мм
r=6.0;#в мм
rad_depth=1.0;#в мм
ax_depth=6.0;#в мм
n=2;
#w_vr=2.19/0.43;#в тысячах оборотов в минуту
#w_vr=2.19/1.75;
#w_vr=2.19/0.63;
w_vr=2.19/1.38;
t_ob=500;#количество шагов по времени на 1 оборот фрезы
eps_F=1.0;#погрешность силы для итерационного алгоритма
Sx=100.0;#в микрометрах
dx_mkm=2.0;#в мкм
dyn_trig=1;#1-динамика включена,0-нет
flank_angle=84;#угол задней грани в градусах
fl_w=1.0;#ширина задней грани в мм
eps_h=0.53;#отношение допустимого проникновения задней грани в конце уточнения к начальному


#пересчет величин
yz=-rad_depth;#в мм
x0=r;
y0=r+yz;
w=1000*w_vr*pi/30;#в рад/с
V=Sx/1000*n*w/(2*pi);
dx=dx_mkm/1000;

#подсчет величин для задней грани
flank_angle_rad=flank_angle*pi/180;#перевод в радианы
r_fl=sqrt(r**2+fl_w**2-2*r*fl_w*cos(flank_angle_rad));
fl_ang=arcsin(fl_w/r_fl*sin(flank_angle_rad));

#коэффициенты сил резания
Krc=445.0;
Ktc=892.0;

#коэффициенты для определения сил, действующих по задней грани
Kn_fl0=1.0*Krc;
mu_fl=0.3;
#Kt_fl=0.0;

#ширина резца
b_rez=ax_depth;#осевая глубина резания

#модальные характеристики
#K_x=2.3*10**6;
K_x=2.3*10**6;
p_x=73*2*pi;
n_x=0.0057*p_x;
m_x=K_x/(p_x**2);
p1_x=sqrt(p_x**2-n_x**2);

#K_y=2.0*10**8;
K_y=2.0*10**8;
p_y=682*2*pi;
n_y=0.0057*p_y;
m_y=K_y/(p_y**2);
p1_y=sqrt(p_y**2-n_y**2);

#определение времени и шагов
step=int(L/dx);
#time=float(L-4*r)/V;

w_num=400;
time_step=int(w_num*t_ob);
dt=2*pi/(w*t_ob);

#time=10.0;
#dt=2*pi/(w*t_ob);
#time_step=int(time/dt);

#массив координат поверхности
z=zeros(step);

#массив координат вершин резцов(текущих)
rez=zeros((n,6));
for i in range(n):
    rez[i][0]=i;
    rez[i][1]=2*i*pi/n;
    rez[i][2]=x0-r*sin(rez[i][1]);
    rez[i][3]=y0+r*cos(rez[i][1]);
    rez[i][4]=x0-r_fl*sin(rez[i][1]-fl_ang);
    rez[i][5]=y0+r_fl*cos(rez[i][1]-fl_ang);


#координаты центра фрезы
xc=x0;
yc=y0;

#массив толщины срезаемого слоя для каждой режущей кромки
h=zeros((n,time_step));

#массив окончательных проникновений задней грани для каждой режущей кромки
h_fl=zeros((n,time_step));

#массив окончательных площадей, ометаемых задней гранью каждой режущей кромки
s_fl=zeros((n,time_step));

#массив начальных площадей, ометаемых задней гранью
s_fl0=zeros((time_step));

#массив начальных проникновений задней грани для каждой режущей кромки
h_fl0_i=zeros((n,time_step));
h_fl0=zeros((time_step));

#массив толщины срезаемого слоя для итерационного уточнения
h_fin=zeros((n));

#массив максимальной величины проникновения задней грани для итерационного уточнения
h_fl_fin=zeros((n));

#массив площади, ометаемой задней гранью
s_fl_fin=zeros((n));

#массивы сил в локальных координатах для каждой кромки
Fr=zeros((n));
Ft=zeros((n));

#массивы сил в локальных координатах, действующих по задней грани, для каждой кромки
Fn_fl=zeros((n));
Ft_fl=zeros((n));

#массивы сил в локальных координатах для каждой кромки в конце шага
Fr_fin=zeros((n));
Ft_fin=zeros((n));

#массивы сил в локальных координатах, действующих по задней грани, для каждой кромки в конце шага
Fn_fl_fin=zeros((n));
Ft_fl_fin=zeros((n));

#массивы суммарных сил в глобальных координатах по всем кромкам
Fx=zeros((time_step));
Fy=zeros((time_step));

#массивы сил от резания в глобальных координатах по всем кромкам
Fxc=zeros((time_step));
Fyc=zeros((time_step));

#массивы сил от задней грани в глобальных координатах по всем кромкам
Fxfl=zeros((time_step));
Fyfl=zeros((time_step));

#массив, содержащий координаты вершин резцов
coord=zeros((n,2,time_step));

#массив, содержащий координаты центра фрезы
coord_c=zeros((2,time_step));

#массив-очередь, хранящий положения резцов
data_s=int(pi/(n*w*dt));
rez_pol=zeros((n,2,data_s));
for i in range(data_s):
    for j in range(n):
        rez_pol[j][0][i]=copy(rez[j][2]);
        rez_pol[j][1][i]=copy(rez[j][3]);


#массив-очередь, хранящий положения резцов - для задней грани
rez_flank=zeros((n,3,data_s));
for i in range(n):
    rez_flank[i][0][:]=(rez[i][2]);
    rez_flank[i][1][:]=(-1.05);
    rez_flank[i][2][:]=(0.0);
rez_flank_last_data=zeros((n,3));
#rez_flank=array([0],ndmin=3);
#массив положений резцов, используемый в итерационном уточнении
rez_fin=zeros((n,6));

#массив текущих положений резцов, обрабатываемый на текущем шаге
rez_w=zeros((n,2));

#массив "следующих" положений резцов, обрабатываемый на текущем шаге
rez_w1=zeros((n,2));

#массив динамических перемещений
x_dyn=zeros((time_step));
y_dyn=zeros((time_step));

#массив, хранящий "триггеры" записи в массивы для задней грани по каждой из кромок
fl_arr_trig=zeros((n));

#массив сил от задней грани, хранящихся для итерационного процесса
F_fl_sum=zeros((data_s));



#"искусственное" начало


xc=17500*dx;
x0=xc;
#смещение центра фрезы
xc=xc+V*dt;
yc=yc;
#запись положений центра фрезы
coord_c[0][0]=xc;
coord_c[1][0]=yc;         
for i in range(n):
    # изменение координат решин резцов
    rez[i][1]=rez[i][1]+w*dt;
    rez[i][2]=xc-r*sin(rez[i][1]);
    rez[i][3]=yc+r*cos(rez[i][1]);
    rez[i][4]=xc-r_fl*sin(rez[i][1]-fl_ang);
    rez[i][5]=yc+r_fl*cos(rez[i][1]-fl_ang);
    #запись положений вершин резцов
    coord[i][0][0]=rez[i][2];
    coord[i][1][0]=rez[i][3];
    
    rez_flank[i][0][0]=rez[i][2];
    rez_flank[i][1][0]=-1.05;
    rez_flank[i][2][0]=rez[i][1];
#изменение поверхности
for i in range(20000):
    z[i]=-1.05;


#сила на конец шага
Fx_fin1=0;
Fy_fin1=0;
Fxc_fin=0;
Fyc_fin=0;
Fxfl_fin=0;
Fyfl_fin=0;



#массив для хранения начальных условий
init_cond=zeros((4));#0-x_0,1-Vx_0,2-y_0,3-Vy_0

#массив для хранения данных из функции вычисления динамических перемещений
data_vect_x=zeros((2));
data_vect_y=zeros((2));

#"выключатель"
break_trig=0;

#скорость вращения режущиз кромок
V_w=w*r;

sum_h_fl=0.0;

sum_s_fl=0;

#начало главного цикла
t=1;
trig_fl_loop=0;


sum_s_otn=zeros((time_step));
s_fl_pred=zeros((n,data_s));

K_fl_fin=zeros((time_step));


#51100
#51060
while t<time_step:
    
    curr_data=int(t%data_s);
    curr_data_1=int((t+1)%data_s);
    
    #блок хранение данных на случай итеарционного уточнения по задней грани

    rez_temp=copy(rez);
    h_fin_temp=copy(h_fin);
    s_fl_fin_temp=copy(s_fl_fin);
    h_fl_fin_temp=copy(h_fl_fin);
    rez_pol_temp=copy(rez_pol);
       
    Fx_fin1_temp=copy(Fx_fin1);
    Fy_fin1_temp=copy(Fy_fin1);
                                
    rez_flank_temp=copy(rez_flank);
    
    init_cond_temp=copy(init_cond);
    
    #конец блока
                                         
                                                                                                          
                                                                                                                                                                                                
                                                                                     
                                                                                                                
                                                                                                                                                                      
    for i in range(n):
                  
        #определение толщины срезаемого слоя
        if (rez[i][2]<0) or (abs(x_dyn[t-1])>1.0):
            break_trig=1;
            print "beda prishla";
            break;
        else:
            h[i][t]=h_fin[i];
            h_fl[i][t]=h_fl_fin[i];
            s_fl[i][t]=s_fl_fin[i];
        
        #блок изменения поверхности полоборота назад    
        
        #задание массивов положений резцов текущих и следующих
        #для изменения поверхности
        rez_w[i][0]=copy(rez_pol[i][0][curr_data]);
        rez_w[i][1]=copy(rez_pol[i][1][curr_data]);
        rez_w1[i][0]=copy(rez_pol[i][0][curr_data_1]);
        rez_w1[i][1]=copy(rez_pol[i][1][curr_data_1]);
        
        
        
        
        zi_w=int(rez_w[i][0]/dx);
        zi_w1=int(rez_w1[i][0]/dx);
        
        z1=min(zi_w,zi_w1);
        z2=max(zi_w,zi_w1);
        
        k_r=k_xy(rez_w[i][0],rez_w[i][1],rez_w1[i][0],rez_w1[i][1]);
        b_r=b_xy(rez_w[i][0],rez_w[i][1],rez_w1[i][0],rez_w1[i][1]);
            
        #изменение поверхности
        for j in range(z1+1,z2+1):
            z_pov=y_x(k_r,b_r,j*dx); 
            if (z[j]>z_pov):
                z[j]=z_pov;

        #расчет суммарных сил в глобальных координатах
        Fx[t]=Fx_fin1;
        Fy[t]=Fy_fin1;
        Fxc[t]=Fxc_fin;
        Fyc[t]=Fyc_fin;
        Fxfl[t]=Fxfl_fin;
        Fyfl[t]=Fyfl_fin;
        
    if break_trig==1:
        print "beda!";
        break;
    
    
    ###
    ###блок взаимодействия по задней грани
    ###
    
    
    #проверка условия (1) для записи в массив для задней грани
    for i in range(n):
        if ((rez[i][1]%(2*pi)>=(pi/2)) and (rez[i][1]%(2*pi)<=(3*pi/2))):
            fl_arr_trig[i]=1;
            rez_flank_last_data[i][0]=rez_flank[i][0][data_s-1];
            rez_flank_last_data[i][1]=rez_flank[i][1][data_s-1];
            rez_flank_last_data[i][2]=rez_flank[i][2][data_s-1];
            rez_flank[i][0]=insert(rez_flank[i][0],0,0.0)[0:data_s];
            rez_flank[i][1]=insert(rez_flank[i][1],0,0.0)[0:data_s];
            rez_flank[i][2]=insert(rez_flank[i][2],0,0.0)[0:data_s];

        else:
            fl_arr_trig[i]=0;
            
    
    
    Kn_fl=0.0*Kn_fl0;            
        
        
    
    #итерационное уточнение по силе на конце шага
    iter_trig=0;
    iter_F=0;
    while (iter_trig==0):
        if (iter_F==0):
            #нулевое приближение по силе в конце шага
            Fx_fin0=Fx[t];
            Fy_fin0=Fy[t];    
        elif (iter_F>0):
            #последнее приближение становится текущим
            Fx_fin0=Fx_fin1;
            Fy_fin0=Fy_fin1;
        #поиск динамических перемещений
        data_vect_x=dyn_per(a_F(Fx[t],Fx_fin0,dt),Fx[t],init_cond[0],init_cond[1],n_x,p1_x,m_x,dt);
        data_vect_y=dyn_per(a_F(Fy[t],Fy_fin0,dt),Fy[t],init_cond[2],init_cond[3],n_y,p1_y,m_y,dt);
        if (dyn_trig==0):
            xc_fin=x0+V*dt*(t+1);
            yc_fin=y0;
        elif (dyn_trig==1):
            xc_fin=1000*data_vect_x[0]+x0+V*dt*(t+1);
            yc_fin=1000*data_vect_y[0]+y0;
        Fx_fin1=0;
        Fy_fin1=0;
        Fxc_fin=0;
        Fyc_fin=0;
        Fxfl_fin=0;
        Fyfl_fin=0;
        
        
        #поиск толщины срезаемого слоя на конце шага
        for i in range(n):
            rez_fin[i][1]=rez[i][1]+w*dt;
            rez_fin[i][2]=xc_fin-r*sin(rez_fin[i][1]);
            rez_fin[i][3]=yc_fin+r*cos(rez_fin[i][1]);
            rez_fin[i][4]=xc_fin-r_fl*sin(rez_fin[i][1]-fl_ang);
            rez_fin[i][5]=yc_fin+r_fl*cos(rez_fin[i][1]-fl_ang);
            h_fin[i]=h_det(z,dx,rez_fin[i][2],rez_fin[i][3],xc_fin,yc_fin);
            Fr_fin[i]=Krc*h_fin[i]*b_rez;
            Ft_fin[i]=Ktc*h_fin[i]*b_rez;
            
            Fx_fin1=Fx_fin1+Fr_fin[i]*sin(rez_fin[i][1])+Ft_fin[i]*cos(rez_fin[i][1]);
            Fy_fin1=Fy_fin1-Fr_fin[i]*cos(rez_fin[i][1])+Ft_fin[i]*sin(rez_fin[i][1]);

            #поиск максимального проникновения задней грани в поверхность
            if (fl_arr_trig[i]==1):#проверка условия (1)
                h_fl_fin[i]=h_fl_det(rez_fin[i],rez_flank[i],data_s,xc_fin,yc_fin);
                s_fl_fin[i]=s_fl_det(rez_fin[i],rez_flank[i],data_s,xc_fin,yc_fin);
            #вычисление сил, действующих по задней грани
            
            ###
            ###
            ###
            if t>0:
                Fn_fl_fin[i]=Kn_fl*s_fl_fin[i]*b_rez;
                Ft_fl_fin[i]=mu_fl*Fn_fl_fin[i];
            
            
            ###
            ###
            ###
            
            ###добавление сил от задней грани
            #Fx_fin1=Fx_fin1-Fn_fl_fin[i]*cos(rez_fin[i][1]+flank_angle_rad);
            #Fy_fin1=Fy_fin1-Fn_fl_fin[i]*sin(rez_fin[i][1]+flank_angle_rad);
            Fxc_fin=Fxc_fin+Fr_fin[i]*sin(rez_fin[i][1])+Ft_fin[i]*cos(rez_fin[i][1]);
            Fyc_fin=Fyc_fin-Fr_fin[i]*cos(rez_fin[i][1])+Ft_fin[i]*sin(rez_fin[i][1]);
            Fxfl_fin=Fxfl_fin-Fn_fl_fin[i]*cos(rez_fin[i][1]+flank_angle_rad)+Ft_fl_fin[i]*sin(rez_fin[i][1]+flank_angle_rad);
            Fyfl_fin=Fyfl_fin-Fn_fl_fin[i]*sin(rez_fin[i][1]+flank_angle_rad)-Ft_fl_fin[i]*cos(rez_fin[i][1]+flank_angle_rad);
            
            Fx_fin1=Fx_fin1-Fn_fl_fin[i]*cos(rez_fin[i][1]+flank_angle_rad)+Ft_fl_fin[i]*sin(rez_fin[i][1]+flank_angle_rad);
            Fy_fin1=Fy_fin1-Fn_fl_fin[i]*sin(rez_fin[i][1]+flank_angle_rad)-Ft_fl_fin[i]*cos(rez_fin[i][1]+flank_angle_rad);


        iter_F=iter_F+1;
        if (iter_F>100):
            print "iter_F >100 1";
            break;
        #сравнение двух "соседних" приближений
        if (sqrt((Fx_fin1-Fx_fin0)**2+(Fy_fin1-Fy_fin0)**2)<=eps_F):
            iter_trig=1;
        
        ###
        ###
        ###
        
    
    
    
    
    ###
    ###блок задней грани
    ###
    
    
    ###формирование поверхности
    for i in range(n):
        #определение скорости вершины режущей кромки
        V_dyn_x=1000*data_vect_x[1];
        V_dyn_y=1000*data_vect_y[1];
        V_rez=V_w-V_dyn_x*cos(rez[i][1])-V_dyn_y*sin(rez[i][1]);
        #проверка на выполнение условия (2) для записи в массив для задней грани
        if (fl_arr_trig[i]==1):#проверка условия (1)
            if (V_rez>=0):#проверка условия (2) - "да"
                ###определение значения, которое будет записываться в массив
                if (h_fin[i]>0):
                    fl_data0=rez_fin[i][2];
                    fl_data1=rez_fin[i][3];
                    fl_data2=rez_fin[i][1];
                elif (h_fin[i]==0):
                    z_curr=int(rez_fin[i][2]/dx);
                    k_curr_s=k_xy(z_curr*dx,z[z_curr],(z_curr+1)*dx,z[z_curr+1]);
                    b_curr_s=b_xy(z_curr*dx,z[z_curr],(z_curr+1)*dx,z[z_curr+1]);
                    fl_data0=rez_fin[i][2];
                    fl_data1=y_x(k_curr_s,b_curr_s,rez_fin[i][2]);
                    fl_data2=rez_fin[i][1];
                if (rez_flank[i][0][1]<=rez_fin[i][2]) or (int(rez_fin[i][1]/(2*pi))>=int(rez_flank[i][2][data_s-1]/(2*pi))):#проверка условия (3)- "да" 
                    ###запись в массив соответствующей точки(положение вершины резца, если режет,
                    ###иначе - соответствующая точка на поверхности
                    rez_flank[i][0][0]=fl_data0;
                    rez_flank[i][1][0]=fl_data1;
                    rez_flank[i][2][0]=fl_data2;
                else:#"нет" для условия (3)
                    #запись
                    rez_flank[i][0][0]=fl_data0;
                    rez_flank[i][1][0]=fl_data1;
                    rez_flank[i][2][0]=fl_data2;
                    #сортировка
                    sorted_index=argsort(-rez_flank[i][0]);
                    rez_flank[i][0]=rez_flank[i][0][sorted_index];
                    rez_flank[i][1]=rez_flank[i][1][sorted_index];
                    rez_flank[i][2]=rez_flank[i][2][sorted_index];         
            else:#проверка условия (2) - "нет"
                rez_flank[i][0]=insert(rez_flank[i][0][1:data_s],data_s-1,rez_flank_last_data[i][0]);
                rez_flank[i][1]=insert(rez_flank[i][1][1:data_s],data_s-1,rez_flank_last_data[i][1]);

    
        
    
        
        
    
    
    
    #####################
    #####################
    #####################
    #####################
    #####################
    
    
        
    if max(s_fl_fin)!=0:
        t_start=t;
        trig_fl_loop=1;
        iter_fl=0;
        t_fin=time_step;
        
        sum_s_fl0=0;
        
        sliding_trig=0;
        
        F_fl_sum[:]=0;
        for i in range(n):
            s_fl_pred[i][:]=0;
        
        
    else:
        if (t%10000==0):
            print t;
        x_dyn[t]=1000*data_vect_x[0];
        y_dyn[t]=1000*data_vect_y[0];
        #смещение центра фрезы
        
        if (dyn_trig==0):
            xc=x0+V*dt*(t+1);
            yc=y0;
        elif (dyn_trig==1):
            xc=x_dyn[t]+x0+V*dt*(t+1);
            yc=y0+y_dyn[t];
        
            
        init_cond[0]=data_vect_x[0];
        init_cond[1]=data_vect_x[1];
        init_cond[2]=data_vect_y[0];
        init_cond[3]=data_vect_y[1];
        
        #запись положений центра фрезы
        coord_c[0][t]=xc;
        coord_c[1][t]=yc;         
        
        for i in range(n):
            
            #заполнение текущей ячейки массива "пол-оборота назад"                                                        
            rez_pol[i][0][curr_data]=copy(rez[i][2]);
            rez_pol[i][1][curr_data]=copy(rez[i][3]);
            
            
            # изменение координат решин резцов
            rez[i][1]=rez[i][1]+w*dt;
            rez[i][2]=xc-r*sin(rez[i][1]);
            rez[i][3]=yc+r*cos(rez[i][1]);
            rez[i][4]=xc-r_fl*sin(rez[i][1]-fl_ang);
            rez[i][5]=yc+r_fl*cos(rez[i][1]-fl_ang);
            #запись положений вершин резцов
            coord[i][0][t]=rez[i][2];
            coord[i][1][t]=rez[i][3];
            
        t=t+1;
        
                
        
    while trig_fl_loop==1:
        
        #проверка условия "проскальзывания"
        if sliding_trig==1:
            trig_fl_loop=0;
            print "perebor!",iter_fl,t_start,t_fin;
        
        
        #блок инициации начала итерационного цикла
        
        
        rez=copy(rez_temp);
        h_fin=copy(h_fin_temp);
        s_fl_fin=copy(s_fl_fin_temp);
        h_fl_fin=copy(h_fl_fin_temp);
        rez_pol=copy(rez_pol_temp);
        Fx_fin1=copy(Fx_fin1_temp);
        Fy_fin1=copy(Fy_fin1_temp);
        rez_flank=copy(rez_flank_temp);
        init_cond=copy(init_cond_temp);
        
        
        
        #конец блока
        
                
        t=t_start;
        sum_s_fl=0;
        
        
        
        ######
        ######
        ######
        
        while t<=t_fin:
                   
            curr_data=int(t%data_s);
            curr_data_1=int((t+1)%data_s);                       
                                    
            for i in range(n):
                        
                #определение толщины срезаемого слоя
                if rez[i][2]<0:
                    break_trig=1;
                    print "beda prishla1";
                    break;
                else:
                    h[i][t]=h_fin[i];
                    if iter_fl==0:
                        s_fl0[t]=max(s_fl_fin);
                    s_fl[i][t]=s_fl_fin[i];
                    h_fl[i][t]=h_fl_fin[i];
                
                #блок изменения поверхности полоборота назад    
                
                #задание массивов положений резцов текущих и следующих
                #для изменения поверхности
                rez_w[i][0]=copy(rez_pol[i][0][curr_data]);
                rez_w[i][1]=copy(rez_pol[i][1][curr_data]);
                rez_w1[i][0]=copy(rez_pol[i][0][curr_data_1]);
                rez_w1[i][1]=copy(rez_pol[i][1][curr_data_1]);
                
                
                
                
                zi_w=int(rez_w[i][0]/dx);
                zi_w1=int(rez_w1[i][0]/dx);
                
                z1=min(zi_w,zi_w1);
                z2=max(zi_w,zi_w1);
                
                k_r=k_xy(rez_w[i][0],rez_w[i][1],rez_w1[i][0],rez_w1[i][1]);
                b_r=b_xy(rez_w[i][0],rez_w[i][1],rez_w1[i][0],rez_w1[i][1]);
                    
                #изменение поверхности
                for j in range(z1+1,z2+1):
                    z_pov=y_x(k_r,b_r,j*dx); 
                    if (z[j]>z_pov):
                        z[j]=z_pov;
        
                #расчет суммарных сил в глобальных координатах
                if t!=t_start:
                    Fx[t]=Fx_fin1;
                    Fy[t]=Fy_fin1;
                    Fxc[t]=Fxc_fin;
                    Fyc[t]=Fyc_fin;
                    Fxfl[t]=Fxfl_fin;
                    Fyfl[t]=Fyfl_fin;
                               
            if break_trig==1:
                print "beda!";
                break;
            
            
            ###
            ###блок взаимодействия по задней грани
            ###
            
            
            #проверка условия (1) для записи в массив для задней грани
            for i in range(n):
                if ((rez[i][1]%(2*pi)>=(pi/2)) and (rez[i][1]%(2*pi)<=(3*pi/2))):
                    fl_arr_trig[i]=1;
                    rez_flank_last_data[i][0]=rez_flank[i][0][data_s-1];
                    rez_flank_last_data[i][1]=rez_flank[i][1][data_s-1];
                    rez_flank_last_data[i][2]=rez_flank[i][2][data_s-1];
                    rez_flank[i][0]=insert(rez_flank[i][0],0,0.0)[0:data_s];
                    rez_flank[i][1]=insert(rez_flank[i][1],0,0.0)[0:data_s];
                    rez_flank[i][2]=insert(rez_flank[i][2],0,0.0)[0:data_s];
        
                else:
                    fl_arr_trig[i]=0;
                    
            
                  
            
            Kn_fl=0.0;
            
                      
            
            #итерационное уточнение по силе на конце шага
            iter_trig=0;
            iter_F=0;
            while (iter_trig==0):
                if (iter_F==0):
                    #нулевое приближение по силе в конце шага
                    Fx_fin0=Fx[t];
                    Fy_fin0=Fy[t];    
                elif (iter_F>0):
                    #последнее приближение становится текущим
                    Fx_fin0=Fx_fin1;
                    Fy_fin0=Fy_fin1;
                #поиск динамических перемещений
                data_vect_x=dyn_per(a_F(Fx[t],Fx_fin0,dt),Fx[t],init_cond[0],init_cond[1],n_x,p1_x,m_x,dt);
                data_vect_y=dyn_per(a_F(Fy[t],Fy_fin0,dt),Fy[t],init_cond[2],init_cond[3],n_y,p1_y,m_y,dt);
                if (dyn_trig==0):
                    xc_fin=x0+V*dt*(t+1);
                    yc_fin=y0;
                elif (dyn_trig==1):
                    xc_fin=1000*data_vect_x[0]+x0+V*dt*(t+1);
                    yc_fin=1000*data_vect_y[0]+y0;
                Fx_fin1=0;
                Fy_fin1=0;
                Fxc_fin=0;
                Fyc_fin=0;
                Fxfl_fin=0;
                Fyfl_fin=0;
                #поиск толщины срезаемого слоя на конце шага
                for i in range(n):
                    rez_fin[i][1]=rez[i][1]+w*dt;
                    rez_fin[i][2]=xc_fin-r*sin(rez_fin[i][1]);
                    rez_fin[i][3]=yc_fin+r*cos(rez_fin[i][1]);
                    rez_fin[i][4]=xc_fin-r_fl*sin(rez_fin[i][1]-fl_ang);
                    rez_fin[i][5]=yc_fin+r_fl*cos(rez_fin[i][1]-fl_ang);
                    h_fin[i]=h_det(z,dx,rez_fin[i][2],rez_fin[i][3],xc_fin,yc_fin);
                    Fr_fin[i]=Krc*h_fin[i]*b_rez;
                    Ft_fin[i]=Ktc*h_fin[i]*b_rez;
                        
                                        
                    Fx_fin1=Fx_fin1+Fr_fin[i]*sin(rez_fin[i][1])+Ft_fin[i]*cos(rez_fin[i][1]);
                    Fy_fin1=Fy_fin1-Fr_fin[i]*cos(rez_fin[i][1])+Ft_fin[i]*sin(rez_fin[i][1]);
                    
                    Fxc_fin=Fxc_fin+Fr_fin[i]*sin(rez_fin[i][1])+Ft_fin[i]*cos(rez_fin[i][1]);
                    Fyc_fin=Fyc_fin-Fr_fin[i]*cos(rez_fin[i][1])+Ft_fin[i]*sin(rez_fin[i][1]);
        
                    #поиск максимального проникновения задней грани в поверхность
                    if (fl_arr_trig[i]==1):#проверка условия (1)
                        h_fl_fin[i]=h_fl_det(rez_fin[i],rez_flank[i],data_s,xc_fin,yc_fin);
                        s_fl_fin[i]=s_fl_det(rez_fin[i],rez_flank[i],data_s,xc_fin,yc_fin);
                    
                    #if s_fl_fin[i]!=0:
                    #    print iter_fl,iter_F,t,s_fl_fin[i];
                    
                    
                    #вычисление сил, действующих по задней грани
                    
                    ###
                    ###
                    ###
                    
                    #if s_fl_fin[i]!=0:
                    #    print iter_fl,iter_F,t%t_start,s_fl_fin[i],s_fl_pred[i][t%t_start];
                        
                        
                    if iter_fl!=0:                        
                        if s_fl_fin[i]!=0:
                            #Kn_fl=K_fl_cont(s_fl_fin[i],Krc);
                            #Kn_fl=K_fl_cont(abs(s_fl_fin[i]-s_fl_pred[i][t%t_start]),Krc);
                            #Kn_fl=K_fl_cont(s_fl_fin[i]-s_fl_pred[i][t%t_start]+1.0*(sum(s_fl_pred[i][(t+1)%t_start:])),Krc);
                            #Kn_fl=K_fl_cont(sum(s_fl_pred[i][t%t_start:]),Krc);
                            
                            #Kn_fl=K_fl_cont(sum(s_fl_pred[i][t%t_start:t_fin%t_start]*(1-(linspace(0,t_fin%t,t_fin-t))/(t_fin-t))),Krc);#зависимость от суммы последующих площадей проникновения с "коэффициентами влияния"
                            Kn_fl=K_fl_cont1(sum(s_fl_pred[i][t%t_start:t_fin%t_start]*(1-(linspace(0,t_fin%t,t_fin-t))/(t_fin-t))),2*sum(s_fl0[t:t_fin]*(1-(linspace(0,t_fin%t,t_fin-t))/(t_fin-t))),Krc);
                            
                            #if s_fl_fin[i]!=0:
                            #    print iter_fl,t%t_start,Kn_fl,s_fl_fin[i],s_fl_pred[i][t%t_start];
                                                    
                            Fn_fl_fin[i]=Kn_fl*s_fl_fin[i]*b_rez;
                            
                            if sliding_trig==0:
                                if s_fl_fin[i]!=0:
                                    Fn_fl_fin[i]=Fn_fl_fin[i]+F_fl_sum[t%t_start];
                            else:
                                Fn_fl_fin[i]=0*K_fl_cont(sum(s_fl0[t:t_fin]*(1-(linspace(0,t_fin%t,t_fin-t))/(t_fin-t))),Krc)*s_fl0[t]*b_rez;
                                
                            Ft_fl_fin[i]=mu_fl*Fn_fl_fin[i];
                            
                        else:
                            Fn_fl_fin[i]=0;
                            Ft_fl_fin[i]=0;
                        
                    
                                        
                    
                    ###
                    ###
                    ###
                    
                    ###добавление сил от задней грани
                    #Fx_fin1=Fx_fin1-Fn_fl_fin[i]*cos(rez_fin[i][1]+flank_angle_rad);
                    #Fy_fin1=Fy_fin1-Fn_fl_fin[i]*sin(rez_fin[i][1]+flank_angle_rad);
                    
                    Fxfl_fin=Fxfl_fin-Fn_fl_fin[i]*cos(rez_fin[i][1]+flank_angle_rad)+Ft_fl_fin[i]*sin(rez_fin[i][1]+flank_angle_rad);
                    Fyfl_fin=Fyfl_fin-Fn_fl_fin[i]*sin(rez_fin[i][1]+flank_angle_rad)-Ft_fl_fin[i]*cos(rez_fin[i][1]+flank_angle_rad);
                    
                    
                    Fx_fin1=Fx_fin1-Fn_fl_fin[i]*cos(rez_fin[i][1]+flank_angle_rad)+Ft_fl_fin[i]*sin(rez_fin[i][1]+flank_angle_rad);
                    Fy_fin1=Fy_fin1-Fn_fl_fin[i]*sin(rez_fin[i][1]+flank_angle_rad)-Ft_fl_fin[i]*cos(rez_fin[i][1]+flank_angle_rad);
        
        
                iter_F=iter_F+1;
                if (iter_F>100):
                    #print "iter_F >100",t;
                    break;
                #сравнение двух "соседних" приближений
                if (sqrt((Fx_fin1-Fx_fin0)**2+(Fy_fin1-Fy_fin0)**2)<=eps_F):
                    iter_trig=1;
                
                ###
                ###
                ###
                
            #"обновление" текущего F_fl_sum для шага t%t_start    
            if iter_fl!=0:
                if max(s_fl_fin)!=0:
                    #print t,t%t_start,Kn_fl*max(s_fl_fin)*b_rez;
                    F_fl_sum[t%t_start]=max(Fn_fl_fin);
                    #print iter_fl,t%t_start,F_fl_sum[t%t_start];
            
            #print Fn_fl_fin,s_fl_fin,K_fl_cont(s_fl_fin[1],Krc),iter_fl,t%t_start;     
            
            if iter_fl==0:
                sum_s_fl0=sum_s_fl0+max(s_fl_fin);
            sum_s_fl=sum_s_fl+max(s_fl_fin);
            
            if iter_fl!=0:
                for i in range(n):
                    s_fl_pred[i][t%t_start]=copy(s_fl_fin[i]);
            
            if iter_fl!=0:
                K_fl_fin[t]=Kn_fl/Krc;
            
            ###
            ###блок задней грани
            ###
            
            
            ###формирование поверхности
            for i in range(n):
                #определение скорости вершины режущей кромки
                V_dyn_x=1000*data_vect_x[1];
                V_dyn_y=1000*data_vect_y[1];
                V_rez=V_w-V_dyn_x*cos(rez[i][1])-V_dyn_y*sin(rez[i][1]);
                #проверка на выполнение условия (2) для записи в массив для задней грани
                if (fl_arr_trig[i]==1):#проверка условия (1)
                    if (V_rez>=0):#проверка условия (2) - "да"
                        ###определение значения, которое будет записываться в массив
                        if (h_fin[i]>0):
                            fl_data0=rez_fin[i][2];
                            fl_data1=rez_fin[i][3];
                            fl_data2=rez_fin[i][1];
                        elif (h_fin[i]==0):
                            z_curr=int(rez_fin[i][2]/dx);
                            k_curr_s=k_xy(z_curr*dx,z[z_curr],(z_curr+1)*dx,z[z_curr+1]);
                            b_curr_s=b_xy(z_curr*dx,z[z_curr],(z_curr+1)*dx,z[z_curr+1]);
                            fl_data0=rez_fin[i][2];
                            fl_data1=y_x(k_curr_s,b_curr_s,rez_fin[i][2]);
                            fl_data2=rez_fin[i][1];
                        if (rez_flank[i][0][1]<=rez_fin[i][2]) or (int(rez_fin[i][1]/(2*pi))>=int(rez_flank[i][2][data_s-1]/(2*pi))):#проверка условия (3)- "да" 
                            ###запись в массив соответствующей точки(положение вершины резца, если режет,
                            ###иначе - соответствующая точка на поверхности
                            rez_flank[i][0][0]=fl_data0;
                            rez_flank[i][1][0]=fl_data1;
                            rez_flank[i][2][0]=fl_data2;
                        else:#"нет" для условия (3)
                            #запись
                            rez_flank[i][0][0]=fl_data0;
                            rez_flank[i][1][0]=fl_data1;
                            rez_flank[i][2][0]=fl_data2;
                            #сортировка
                            sorted_index=argsort(-rez_flank[i][0]);
                            rez_flank[i][0]=rez_flank[i][0][sorted_index];
                            rez_flank[i][1]=rez_flank[i][1][sorted_index];
                            rez_flank[i][2]=rez_flank[i][2][sorted_index];         
                    else:#проверка условия (2) - "нет"
                        print "negative V";
                        rez_flank[i][0]=insert(rez_flank[i][0][1:data_s],data_s-1,rez_flank_last_data[i][0]);
                        rez_flank[i][1]=insert(rez_flank[i][1][1:data_s],data_s-1,rez_flank_last_data[i][1]);
        
            
                
            if (t%10000==0):
                print t;
            x_dyn[t]=1000*data_vect_x[0];
            y_dyn[t]=1000*data_vect_y[0];
            #смещение центра фрезы
            
            if (dyn_trig==0):
                xc=x0+V*dt*(t+1);
                yc=y0;
            elif (dyn_trig==1):
                xc=x_dyn[t]+x0+V*dt*(t+1);
                yc=y0+y_dyn[t];
            
                
            init_cond[0]=data_vect_x[0];
            init_cond[1]=data_vect_x[1];
            init_cond[2]=data_vect_y[0];
            init_cond[3]=data_vect_y[1];
            
            #запись положений центра фрезы
            coord_c[0][t]=xc;
            coord_c[1][t]=yc;         
            
            for i in range(n):
                
                #заполнение текущей ячейки массива "пол-оборота назад"                                                        
                rez_pol[i][0][curr_data]=copy(rez[i][2]);
                rez_pol[i][1][curr_data]=copy(rez[i][3]);
                
                
                # изменение координат решин резцов
                rez[i][1]=rez[i][1]+w*dt;
                rez[i][2]=xc-r*sin(rez[i][1]);
                rez[i][3]=yc+r*cos(rez[i][1]);
                rez[i][4]=xc-r_fl*sin(rez[i][1]-fl_ang);
                rez[i][5]=yc+r_fl*cos(rez[i][1]-fl_ang);
                #запись положений вершин резцов
                coord[i][0][t]=rez[i][2];
                coord[i][1][t]=rez[i][3];
            
            if iter_fl==0:
                if max(s_fl_fin)==0:
                    t_fin=t;
                
            t=t+1;
        
        
        
        
        
        
        
        
        
        
        #print F_fl_sum;
        
        ######
        ######
        ######
        
        iter_fl=iter_fl+1;
        
        #проверка условия "скольжения"
        if (t_fin-t_start)<=10 or iter_fl==999:
            sliding_trig=1;
            #trig_fl_loop=0;
        
 
        #для расчет без учета задней грани    
        #if iter_fl>0:
        #    print t_start,t_fin,sum_s_fl/sum_s_fl0;
        #    trig_fl_loop=0;
        
        if sum_s_fl/sum_s_fl0<eps_h  or iter_fl>1000:
            print t_start,t_fin,sum_s_fl/sum_s_fl0,iter_F,iter_fl;
            sum_s_otn[t_start:t_fin]=sum_s_fl/sum_s_fl0;
            trig_fl_loop=0;
            
            
    
    
    
    ################
    ################
    ################
    ################
    ################
    
    
        

#print sum_h_fl/sum_h_fl0;
#print sum_s_fl/sum_s_fl0;






"""

fig=plt.figure(200);
ax=fig.add_subplot(1,1,1);
ax.set_xlim((0,L));
ax.set_ylim((-2,L-2));


rez_1x,rez_1y=array([[xc,rez_fin[2],rez_fin[4]],[yc,rez_fin[3],rez_fin[5]]]);
line1 = mlines.Line2D(rez_1x,rez_1y);
ax.add_line(line1);
plt.plot(rez_flank[0],rez_flank[1]);

"""

# -*- coding: utf-8 -*-
"""
Created on Sun Nov 27 15:14:22 2016

@author: Xin
"""
#solve 2-D incompressible NS equations using spectral method 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import os.path
from numpy import *
from numpy.fft import fftfreq , fft , ifft , irfft2 , rfft2, fft2, ifft2, fftshift, ifftshift
from mpi4py import MPI
parent = os.path.abspath(os.path.join(os.path.dirname(__file__),'.'))
sys.path.append(parent)


#parameteres
new=1;
forcing=1;

Nstep=5001; #no. of steps
Nx=64; #grid size
Ny=64;
t=0;

nu=5e-10; #viscosity
nu_hypo=2e-3; #hypo-viscosity
dt=5e-4; #time-step
dt_h=dt/2; #half-time step
init_cond = 'IC_TaylorGreen';
k_ic=1;         #for Taylor-Green init_cond
force_cond = 'F_TaylorGreen';
k_f=2;          #for Taylor-Green forcing

diag_out_step = 5000; #step frequency of outputting diagnostics


#-----------GRID setup-----------
Lx=2*pi;
Ly=2*pi;
dx=Lx/Nx;
dy=Ly/Ny;
x=zeros((Nx, Ny), dtype=float);
y=zeros((Nx, Ny), dtype=float);
kx=zeros((Nx, Ny), dtype=float);
ky=zeros((Nx, Ny), dtype=float);
for j in range(Ny):
    kx[0:Nx,j]=range(-Nx/2, Nx/2);
    x[0:Nx,j] =range(-Nx/2, Nx/2);
x=x*dx;
for i in range(Nx):
    ky[i,0:Ny]=range(-Ny/2, Ny/2);
    y[i,0:Ny] =range(-Ny/2, Ny/2);
y=y*dy;
k2=kx**2+ky**2;
for i in range(Nx):
    for j in range(Ny):
         if(k2[i,j] == 0):
             k2[i,j]=1e-5; #so that I do not divide by 0 below when using projection operator

k2_exp=exp(-nu*(k2**5)*dt-nu_hypo*dt);

#-----------Force setup----------
Fxhat = zeros((Nx, Ny), dtype=complex);
Fyhat = zeros((Nx, Ny), dtype=complex);
for iss in [-1, 1]:
    for jss in [-1, 1]:
        for i in range(Nx):
            for j in range(Ny):
                if(int(kx[i,j])==iss*k_f and int(ky[i,j])==jss*k_f):
                    Fxhat[i,j] = -1j*iss;
                    Fyhat[i,j] = -1j*(-jss);
                
Fxhat=5e-2*Fxhat;     
Fyhat=5e-2*Fyhat;

#assume no foricng
Fxhat=0.0*Fxhat
Fyhat=0.0*Fyhat


#-----------Variables------------
Vxhat=zeros((Nx, Ny), dtype=complex);
Vyhat=zeros((Nx, Ny), dtype=complex);
Wzhat=zeros((Nx, Ny), dtype=complex);
Vx=zeros((Nx, Ny), dtype=float);
Vy=zeros((Nx, Ny), dtype=float);
Wz=zeros((Nx, Ny), dtype=float);


#----Initialize: Taylor-Green Velocity in Fourier space-----------
if (new==1):
    Vxhat = zeros((Nx, Ny), dtype=complex);
    Vyhat = zeros((Nx, Ny), dtype=complex);
    for iss in [-1, 1]:
        for jss in [-1, 1]:
            for i in range(Nx):
                for j in range(Ny):
                    if(int(kx[i,j])==iss*k_ic and int(ky[i,j])==jss*k_ic):
                        Vxhat[i,j] = -1j*iss;
                        Vyhat[i,j] = -1j*(-jss);
                
#Set total energy to 1
    Vxhat=0.5*Vxhat;     
    Vyhat=0.5*Vyhat;
    
    Vx=10*random.rand(Nx,Ny)
    Vy=10*random.rand(Nx,Ny)
    
    Vxhat = (fftshift(fft2(ifftshift(Vx))));
    Vyhat = (fftshift(fft2(ifftshift(Vy))));
    
#----read from data---------
else:
    Vxhat = Vxhat_restart;
    Vyhat = Vyhat_restart;


#------Dealiasing------------------------------------------------
for i in range(Nx):
      for j in range(Ny):
          if(sqrt(k2[i,j]) >= Nx/3.):
              #print k2[i,j]
              Vxhat[i,j]=0;
              Vyhat[i,j]=0;


#------Projection operator on velocity fields to make them solenoidal-----
tmp = (kx*Vxhat + ky*Vyhat)/k2;
Vxhat = Vxhat - kx*tmp;
Vyhat = Vyhat - ky*tmp;

#------Storing variables for later use in time integration--------
Vxhat_t0 = Vxhat;
Vyhat_t0 = Vyhat;

#----Main Loop-----------
for istep in range(Nstep):
    #------Dealiasing------------------------------------------------
    for i in range(Nx):
        for j in range(Ny):
            if(np.sqrt(k2[i,j]) >= Nx/3.):
                Vxhat[i,j]=0;
                Vyhat[i,j]=0;
    #Projection operator on velocity
    tmp = (kx*Vxhat + ky*Vyhat)/k2;
    Vxhat = Vxhat - kx*tmp;
    Vyhat = Vyhat - ky*tmp;
    
    #Calculate Vorticity
    Wzhat = 1j*(kx*Vyhat - ky*Vxhat);
    #fileds in x-space
    Vx=fftshift(ifft2(ifftshift(Vxhat))).real;
    Vy=fftshift(ifft2(ifftshift(Vyhat))).real;
    Wz=fftshift(ifft2(ifftshift(Wzhat))).real;
    
    #Fields in Fourier Space
    Vxhat = (fftshift(fft2(ifftshift(Vx))));
    Vyhat = (fftshift(fft2(ifftshift(Vy))));
    Wzhat = (fftshift(fft2(ifftshift(Wz))));
   
    #Calculate non-linear term in x-space
    NLx =  Vy*Wz;
    NLy = -Vx*Wz;
    
    
    #move non-linear term back to Fourier k-space
    NLxhat = (fftshift(fft2(ifftshift(NLx))));
    NLyhat = (fftshift(fft2(ifftshift(NLy))));
 

    
    #------Dealiasing------------------------------------------------
    for i in range(Nx):
        for j in range(Ny):
            if(sqrt(k2[i,j]) >= Nx/3.):
                NLxhat[i,j]=0;
                NLyhat[i,j]=0;
                
    #------Projection operator on velocity fields to make them solenoidal-----
    tmp = (kx*NLxhat + ky*NLyhat)/k2;
    NLxhat = NLxhat - kx*tmp;
    NLyhat = NLyhat - ky*tmp;
    
    #Integrate in time
    #---Euler for 1/2-step-----------
    if(istep==0):
          Vxhat = Vxhat + dt_h*(NLxhat -nu*(k2**5)*Vxhat -nu_hypo*(k2**(-0))**Vxhat) + dt_h**Fxhat;
          Vyhat = Vyhat + dt_h*(NLyhat -nu*(k2**5)*Vyhat -nu_hypo*(k2**(-0))**Vyhat) + dt_h**Fyhat;
          oNLxhat = NLxhat;
          oNLyhat = NLyhat;
    #---Midpoint time-integration----
    elif(istep==1):
          Vxhat = Vxhat_t0 + dt*(NLxhat -nu*(k2**5)*Vxhat -nu_hypo*(k2**(-0))*Vxhat) + dt*Fxhat;
          Vyhat = Vyhat_t0 + dt*(NLyhat -nu*(k2**5)*Vyhat -nu_hypo*(k2**(-0))*Vyhat) + dt*Fyhat;          
    #---Adam-Bashforth integration---
    else:
          Vxhat = Vxhat + dt*(1.5*NLxhat - 0.5*oNLxhat*k2_exp);
          Vyhat = Vyhat + dt*(1.5*NLyhat - 0.5*oNLyhat*k2_exp);
          Vxhat = Vxhat*k2_exp;
          Vyhat = Vyhat*k2_exp;
          Vxhat = Vxhat + dt*Fxhat;
          Vyhat = Vyhat + dt*Fyhat;

          oNLxhat=NLxhat;
          oNLyhat=NLyhat;
          
    if(istep%diag_out_step==0):
        plt.contourf(x, y, Wz)
        plt.show()
    t=t+dt;

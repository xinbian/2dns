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
N=Nx=Ny=64;
t=0;

nu=5e-10; #viscosity
nu_hypo=2e-3; #hypo-viscosity
dt=5e-4; #time-step
dt_h=dt/2; #half-time step
k_ic=2;  #1 for Taylor-Green init_cond; 2 for random init_cond
diag_out_step = 5000; #step frequency of outputting diagnostics

#------------MPI setup---------
comm = MPI.COMM_WORLD
num_processes = comm.Get_size()
rank = comm.Get_rank()
#slab decomposition, split in x in physical space, in y in Fourrier space
Np = Nx/num_processes

#---------modifited 2d FFT and IFFT
Uc_hat = empty(( N, Np) , dtype = complex )
Uc_hatT = empty((Np,N) , dtype = complex )
U_mpi = empty (( num_processes , Np , Np ) , dtype = complex )
def ifftn_mpi(fu, u):
    Uc_hat[:]=fftshift(ifft(ifftshift(fu) , axis=0))
    comm.Alltoall([Uc_hat, MPI.DOUBLE_COMPLEX], [U_mpi, MPI.DOUBLE_COMPLEX])
    Uc_hatT[:] = rollaxis(U_mpi, 1).reshape(Uc_hatT.shape)
    u[:]= fftshift(ifft (ifftshift(Uc_hatT), axis =1))
    return u
def fftn_mpi(u, fu):
    Uc_hatT[:] = fftshift(fft(ifftshift(u), axis=1))
    U_mpi[:] = rollaxis(Uc_hatT.reshape(Np, num_processes, Np), 1)
    comm.Alltoall([U_mpi, MPI.DOUBLE_COMPLEX], [fu, MPI.DOUBLE_COMPLEX])
    fu[:]= fftshift(fft(ifftshift(fu), axis=0))
    return fu  

#-----------GRID setup-----------
Lx=2*pi;
Ly=2*pi;
dx=Lx/Nx;
dy=Ly/Ny;
x=zeros((Np, Ny), dtype=float);
y=zeros((Np, Ny), dtype=float);
kx=zeros((Nx, Np), dtype=float);
ky=zeros((Nx, Np), dtype=float);
#initilize x space coordinate
for j in range(Ny):
    x[0:Np,j]=range(Np);
x=x-(num_processes/2-rank)*Np
x=x*dx;

for i in range(Np):
    y[i,0:Ny] =range(-Ny/2, Ny/2);
y=y*dy;

for j in range(Np):
    kx[0:Nx,j]=range(-Nx/2, Nx/2);
for i in range(Nx):
    ky[i,0:Np]=range(Np);
ky=ky-(num_processes/2-rank)*Np

k2=kx**2+ky**2;
for i in range(Nx):
    for j in range(Np):
         if(k2[i,j] == 0):
             k2[i,j]=1e-5; #so that I do not divide by 0 below when using projection operator

k2_exp=exp(-nu*(k2**5)*dt-nu_hypo*dt);

#---------dealiasing function----
def delias(Vxhat, Vyhat):
    for i in range(Nx):
        for j in range(Np):
            if(sqrt(k2[i,j]) >= Nx/3.):
                Vxhat[i,j]=0;
                Vyhat[i,j]=0;
    #------Projection operator on velocity fields to make them solenoidal-----
    tmp = (kx*Vxhat + ky*Vyhat)/k2;
    Vxhat = Vxhat - kx*tmp;
    Vyhat = Vyhat - ky*tmp;
    return Vxhat, Vyhat


#-----------Variables------------
Vxhat=zeros((Nx, Np), dtype=complex);
Vyhat=zeros((Nx, Np), dtype=complex);
Wzhat=zeros((Nx, Np), dtype=complex);
NLxhat=zeros((Nx, Np), dtype=complex);
NLyhat=zeros((Nx, Np), dtype=complex);
Vx=zeros((Np, Ny), dtype=float);
Vy=zeros((Np, Ny), dtype=float);
Wz=zeros((Np, Ny), dtype=float);
Wz_all=zeros((Nx, Ny), dtype=float);


##----Initialize: Taylor-Green Velocity in Fourier space-----------
if (new==1 and k_ic==1):
    Vxhat = zeros((Nx, Np), dtype=complex);
    Vyhat = zeros((Nx, Np), dtype=complex);
    for iss in [-1, 1]:
        for jss in [-1, 1]:
            for i in range(Nx):
                for j in range(Np):
                    if(int(kx[i,j])==iss*k_ic and int(ky[i,j])==jss*k_ic):
                        Vxhat[i,j] = -1j*iss;
                        Vyhat[i,j] = -1j*(-jss);
#                
##Set total energy to 1
    Vxhat=0.5*Vxhat;     
    Vyhat=0.5*Vyhat;
elif (new==1 and k_ic==2):
    Vx=10*random.rand(Np,Ny)
    Vy=10*random.rand(Np,Ny)
    Vxhat=fftn_mpi(Vx, Vxhat)
    Vyhat=fftn_mpi(Vy, Vyhat)
   
#    Vxhat = (fftshift(fft2(ifftshift(Vx))));
#    Vyhat = (fftshift(fft2(ifftshift(Vy))));
#    
##----read from data---------
else:
    Vxhat = Vxhat_restart;
    Vyhat = Vyhat_restart;
#----output vorticity function----
#def output(Wz, rank):
#    for i in range(1,num_processes):
#        if rank==i:
#            comm.send(Wz, dest=0, tag=i) 
#    if rank==0:
#        for i in range(1,num_processes):
#            Wz_all[i*Np,:]=comm.recv(source=i, tag=i)
#if rank == 0:
##   data = {'a': 7, 'b': 3.14}
#   comm.send(data, dest=1, tag=11)
#elif rank == 1:
#   data = comm.recv(source=0, tag=11)
        
#------Dealiasing------------------------------------------------
Vxhat, Vyhat=delias(Vxhat, Vyhat)
#
#------Storing variables for later use in time integration--------
Vxhat_t0 = Vxhat;
Vyhat_t0 = Vyhat;
#
#----Main Loop-----------
for istep in range(Nstep):
    #------Dealiasing------------------------------------------------
    for i in range(Nx):
        for j in range(Np):
            if(sqrt(k2[i,j]) >= Nx/3.):
                Vxhat[i,j]=0;
                Vyhat[i,j]=0;
    #Projection operator on velocity
    tmp = (kx*Vxhat + ky*Vyhat)/k2;
    Vxhat = Vxhat - kx*tmp;
    Vyhat = Vyhat - ky*tmp;
    
    #Calculate Vorticity
    Wzhat = 1j*(kx*Vyhat - ky*Vxhat);
    #fileds in x-space
    #Vx=fftshift(ifft2(ifftshift(Vxhat))).real;
    Vx=ifftn_mpi(Vxhat,Vx)
    #Vy=fftshift(ifft2(ifftshift(Vyhat))).real;
    Vy=ifftn_mpi(Vyhat,Vy)
    #Wz=fftshift(ifft2(ifftshift(Wzhat))).real;
    Wz=ifftn_mpi(Wzhat,Wz)
    
    #Fields in Fourier Space
    #Vxhat = (fftshift(fft2(ifftshift(Vx))));
    #Vyhat = (fftshift(fft2(ifftshift(Vy))));
    #Wzhat = (fftshift(fft2(ifftshift(Wz))));
    Vxhat=fftn_mpi(Vx, Vxhat)
    Vyhat=fftn_mpi(Vy, Vyhat)
    Wzhat=fftn_mpi(Wz, Wzhat)
   
    #Calculate non-linear term in x-space
    NLx =  Vy*Wz;
    NLy = -Vx*Wz;
    
    
    #move non-linear term back to Fourier k-space
    #NLxhat = (fftshift(fft2(ifftshift(NLx))));
    #NLyhat = (fftshift(fft2(ifftshift(NLy))));
    NLxhat=fftn_mpi(NLx, NLxhat)
    NLyhat=fftn_mpi(NLy, NLyhat)
    
    #------Dealiasing------------------------------------------------
    Vxhat, Vyhat=delias(Vxhat, Vyhat)
    
    #Integrate in time
    #---Euler for 1/2-step-----------
    if(istep==0):
          Vxhat = Vxhat + dt_h*(NLxhat -nu*(k2**5)*Vxhat -nu_hypo*(k2**(-0))**Vxhat);
          Vyhat = Vyhat + dt_h*(NLyhat -nu*(k2**5)*Vyhat -nu_hypo*(k2**(-0))**Vyhat);
          oNLxhat = NLxhat;
          oNLyhat = NLyhat;
    #---Midpoint time-integration----
    elif(istep==1):
          Vxhat = Vxhat_t0 + dt*(NLxhat -nu*(k2**5)*Vxhat -nu_hypo*(k2**(-0))*Vxhat);
          Vyhat = Vyhat_t0 + dt*(NLyhat -nu*(k2**5)*Vyhat -nu_hypo*(k2**(-0))*Vyhat);          
    #---Adam-Bashforth integration---
    else:
          Vxhat = Vxhat + dt*(1.5*NLxhat - 0.5*oNLxhat*k2_exp);
          Vyhat = Vyhat + dt*(1.5*NLyhat - 0.5*oNLyhat*k2_exp);
          Vxhat = Vxhat*k2_exp;
          Vyhat = Vyhat*k2_exp;
          Vxhat = Vxhat;
          Vyhat = Vyhat;

          oNLxhat=NLxhat;
          oNLyhat=NLyhat;
          
    if(istep%diag_out_step==0 and rank==1):
        plt.contourf(x, y, Wz)
        plt.show()
        plt.savefig('myfig')
    t=t+dt;

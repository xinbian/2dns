python 2D NS
============
[![image](https://img.shields.io/travis/xinbian/2dns.svg)](https://travis-ci.org/xinbian/2dns)[![image](https://readthedocs.org/projects/2dns-2/badge/?version=latest)](https://2dns-2.readthedocs.io/en/latest/?badge=latest)[![image](https://codecov.io/gh/xinbian/2dns/branch/master/graph/badge.svg)](https://codecov.io/gh/xinbian/2dns)

In the code, the 2D Navier-Stokes (NS) equation is solved.

-   Free software: MIT license
-   Documentation: <https://2dns-2.readthedocs.io>.

Features
--------
1. The 2D incompressible NS equation is solved using [spectral method](https://en.wikipedia.org/wiki/Spectral_method). The boundary condition is periodic. For now, the gird geometry must be square.
2. Explicit time integration method is used, containing Euler method, mid-point method, and Adam-Bashforth method. The time step needs to be decreased when increasing grid resolution.
3. The code is parallelized with MPI. The code only works with cores of power of 2. Slab decomposition is used. In physical space, the array in y direction is decomposed; In Fourier space, the array in kx direction is decomposed.
4. The functions of 2D FFT and IFFT are based on the code of [spectralDNS](https://github.com/spectralDNS/spectralDNS). 3D functions are modified to 2D.

Usage
----------
1. mpi4py need to be installed first.
2. run the code  /python_2d_ns/python_2d_ns.py with command 'mpirun -n 2 python python_2d_ns.py'. The value after n is the number of cores used in the simulation.
3. You need to modify some parameters to customize the simulation.

| N | ic_type                | k_ic              | dt        | nu        | nu_hypo        |
|-----------|------------------------|-------------------|-----------|-----------|----------------|
| grid size | initial condition (IC) | wavenumber for IC | time step | viscosity | hypo-viscoisty |

Parallel
----------
The code is tested on [NESRSC cori supercomputer](http://www.nersc.gov/users/computational-systems/cori/), up to 1024 cores with N=4096.
![alt text](https://pbs.twimg.com/media/CzhTly9WQAEY1up.jpg:large "parallel")

Sample results
-----------
1. IC: Taylor-Green vortex (the vortex geometry will remain the same)
![](https://pbs.twimg.com/media/CzhYsZAXAAIcZwd.jpg:large)
![](https://pbs.twimg.com/media/CzhYtu7XgAA2aB3.jpg:large)
2. IC: random velocity



Credits
----------

This package was created with [Cookiecutter](https://github.com/audreyr/cookiecutter) and the [audreyr/cookiecutter-pypackage](https://github.com/audreyr/cookiecutter-pypackage)
project template.
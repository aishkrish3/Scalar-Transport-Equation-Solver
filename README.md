## Scalar Transport Equation Solver

### Decription:
This solver aims to solve the unsteady, scalar transport euation in 2D. This equations consists of an unsteady term, diffusion term, convective term and source.
Examples of a scalar would be mass or temperature. The latter can be thought of as a heat sink and a source in a 2D mesh. The code implements a uniform cartesian mesh of [0,1]x[1,0] domain with periodic boundary conditions. It comprises of 32 cells in x direction and 32 cells in y direction. 

The project implements various control volume techniques for solving the equation:
 - Second order linear reconstruction to compute the derivative of scalar in the diffusion term - the method by which the cell-centered, averaged values of conserved quantities are interpolated to cell faces, in order to calculate the left- and right-states needed to compute fluxes, in this case the diffusion flux
 - 3rd order QUICK (Quadratic Upstream Interpolation for Convective Kinematics) scheme for convective flux - higher-order differencing scheme that considers a three-point upstream weighted quadratic interpolation for the cell face values
 The project uses two time integration schemes:
    - 1st order explicit
    - 2nd order Crank-Nicolson - It is a second-order implicit method in time which is proven to be unconditionally stable
    
 The time step (Δt) is restricted by stability limit of the solver (i.e., time step is limited by the Courant–Friedrichs–Lewy condition equal to 0.8 in this solver).
    
The plot shows the unsteady transport of the scalar, φ due to the convective, diffusion and source terms. The convective term specifies movement of φ based on velocities. The diffusion term specifies the “spread” of φ.

### How to build:
You will need the following packages in a Debian/Ubuntu environment:
 - GNU C++ compiler with C++17 support
 - Armadillo: http://arma.sourceforge.net/download.html
 - Boost: https://www.boost.org/users/download/
 - GNUPlot: http://www.gnuplot.info/download.html

#### Manual build:
In order to build this project from the command line:
```console
$ cd Debug
$ make
```

#### Automated build:
This repository contains Unix specific Eclipse CDT settings for convenience. 
 
### Author:
Aishwarya Krishnan

## Scalar Transport Equation Solver

### Decription:
This solver aims to solve the unsteady, scalar transport euation in 2D. This equations consists of an unsteady term, diffusion term, convective term and source.
Examples of a scalar would be mass or temperature. The latter can be thought of as a heat sink and a source in a 2D mesh. The code implements a uniform cartesian mesh of [0,1]x[1,0] domain with periodic boundary conditions. It comprises of 32 cells in x direction and 32 cells in y direction. 

The project implements various control volume techniques for solving the equation:
 - Second order linear reconstruction to compute the derivative of scalar in the diffusion term - the method by which the cell-centered, averaged values of conserved quantities are interpolated to cell faces, in order to calculate the left- and right-states needed to compute fluxes, in this case the diffusion flux
 - 3rd order QUICK (Quadratic Upstream Interpolation for Convective Kinematics) scheme for convective flux - higher-order differencing scheme that considers a three-point upstream weighted quadratic interpolation for the cell face values
 - Two time integration schemes:
    - 1st order explicit
    - 2nd order Crank-Nicolson - It is a second-order implicit method in time which is proven to be unconditionally stable
    
 The time step (Δt) is restricted by stability limit of the solver (i.e., time step is limited by the Courant–Friedrichs–Lewy condition equal to 0.8 in this solver).
    
The plot shows the unsteady transport of the scalar, φ due to the convective, diffusion and source terms. The convective term specifies movement of φ based on velocities. The diffusion term specifies the “spread” of φ.

### Installation:
I have used the following dependencies which need to be installed before running the code:
 - Armadillo: http://arma.sourceforge.net/download.html
 - Boost: https://www.boost.org/users/download/
 - GNUPlot: http://www.gnuplot.info/download.html
 
 The code was written using Eclipse on Linux. Once the code has been downloaded, any C++ compiler can be used to run SampleCode.cpp.
 
 ### Class Structure:
 The following is the structure of the classes and their usage:
- SampleCode.cpp contains main() and calls the time integration codes( "CrankNicolson analysis(MESH_SIZE, MESH_SIZE);" for Crank Nicolson and "FirstOrderExplicit analysis(MESH_SIZE, MESH_SIZE);" for the explicit time integration) and Analysis.cpp
- CrankNicolson/FirstOrderExplicit creates the mesh according to MESH_SIZE and calculates convection, diffusion terms
- Analysis.cpp calculates the scalar φ from convection, diffusion terms and solves it using the Armadillo library(this equation is solved differently for both time integration steps based on the techniques specified in the "Description" section)
- Finally, SampleCode.cpp plots the scalar at different timesteps using GNUPlot

### Author:
Aishwarya Krishnan

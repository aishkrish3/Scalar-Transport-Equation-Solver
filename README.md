# Scalar Transport Equation Solver

## Decription:
This project aims to solve an unsteady, scalar transport equation for a uniform, 2D cartesian mesh with periodic boundary conditions. These equations consist of an unsteady term, a diffusion term, a convective term and a source.

The project implements following control volume schemes for solving the equation:
 - [Second order linear reconstruction](https://en.wikipedia.org/wiki/MUSCL_scheme) to compute the derivative of scalar in the diffusion term - the method by which the cell-centered, averaged values of conserved quantities are interpolated to cell faces, in order to calculate the left-state and right-state needed to compute fluxes, in this case the diffusion flux.
 - 3rd order [QUICK](https://en.wikipedia.org/wiki/QUICK_scheme) (Quadratic Upstream Interpolation for Convective Kinematics) scheme for convective flux - higher-order differencing scheme that considers a three-point upstream weighted quadratic interpolation for the cell face values.
 The project uses two [time integration schemes](https://en.wikipedia.org/wiki/Temporal_discretization) specified by the user:
    - 1st order explicit - implemented by the FTCS (Forward Time Central Space) scheme which is conditionally stable.
    - 2nd order Crank-Nicolson - It is a second-order implicit method in time which is proven to be unconditionally stable.
    
 The time step (**_Δt_**) is restricted by stability limit of the solver (i.e., time step is limited by the [Courant–Friedrichs–Lewy number](https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition) equal to 0.8 in this solver).
    
The following plot generated by the program shows the unsteady transport of the scalar, **_φ_** due to the convective, diffusion and source terms using Crank Nicolson at regular timesteps:

![alt_text](https://github.com/aishkrish3/Scalar-Transport-Equation-Solver/blob/master/documentation/Output.png "Scalar transport at regular timesteps using CN time integration")

The convective term specifies movement of **_φ_** based on velocities. The diffusion term specifies the “spread” of **_φ_**. 

## How to build:
You will need the following packages in a Debian/Ubuntu environment:
 - GNU C++ compiler with C++17 support
 - [Armadillo](http://arma.sourceforge.net/download.html), a linear algebra solver for C++.
 - [GNUPlot](http://www.gnuplot.info/download.html), for data visualization.
 - [GNUplot-iostream interface](https://code.google.com/archive/p/gnuplot-cpp/), a C++ interface for GNUPlot.
   - Requires [boost](https://www.boost.org/users/download/) C++ libraries.

### Manual build:
In order to build this project from the command line:
```console
$ cd Debug
$ make
```

### Automated build:
This repository contains Unix specific Eclipse CDT settings for convenience. 

## Running the program:
Building the program generates an executable in the Debug folder.
You can run it using:
```console
$ cd Debug
$ ./ScalarTransportEquation
```

## Theory:
The scalar transport equation is a combination of the diffusion and convection (advection) equations, and describes a physical phenomena where particles, energy, or other physical quantities are transferred inside a physical system due to two processes: diffusion and convection. In its full unsteady form, it is given by:

![alt_text](https://github.com/aishkrish3/Scalar-Transport-Equation-Solver/blob/master/documentation/Scalar%20equation.JPG "Scalar transport Equation")

 - The **diffusion term** is given by **_∇ ⋅ (α∇φ)_**. Imagine that **_φ_** is the concentration of a chemical. When concentration is low somewhere compared to the surrounding areas (e.g. a local minimum of concentration), the substance will diffuse in from the surroundings, so the concentration will increase. Conversely, if concentration is high compared to the surroundings (e.g. a local maximum of concentration), then the substance will diffuse out and the concentration will decrease. The net diffusion is proportional to the Laplacian (or second derivative) of concentration if the diffusivity **_α_** is a constant.
 - The **convective(or advective) term** is **_−∇ ⋅ (vφ)_**. Imagine standing on the bank of a river, measuring the water's salinity (amount of salt) each second. Upstream, somebody dumps a bucket of salt into the river. A while later, you would see the salinity suddenly rise, then fall, as the zone of salty water passes by. Thus, the concentration at a given location can change because of the flow.
 - The **reactive source term** is given by **_ω_** which describes the creation or destruction of the quantity. For example, if **_φ_** is the concentration of a molecule, then **_ω_** describes how the molecule can be created or destroyed by chemical reactions. **_ω_** may be a function of **_φ_** and of other parameters. Often there are several quantities, each with its own convection–diffusion equation, where the destruction of one quantity entails the creation of another. For example, when methane burns, it involves not only the destruction of methane and oxygen but also the creation of carbon dioxide and water vapor. Therefore, while each of these chemicals has its own convection–diffusion equation, they are coupled together and must be solved as a system of simultaneous differential equations.

This project solves for the scalar field **_φ_** in 2D on a [0,1]x[1,0] domain with periodic boundary conditions (i.e., **_φ(x=1,y) = φ(x=0,y) and φ(x,y=1) = φ(x,y=0)_**).

The theoretical **velocity field** is prescribed by:

![alt_text](https://github.com/aishkrish3/Scalar-Transport-Equation-Solver/blob/master/documentation/u_velocity.JPG "Velocity in x direction")

and

![alt_text](https://github.com/aishkrish3/Scalar-Transport-Equation-Solver/blob/master/documentation/v_velocity.JPG "Velocity in y direction")

The **diffusion coefficient** is given by **_α = 1e-2_**.
The theoretical **reaction source (or sink) term** is given by:

![alt_text](https://github.com/aishkrish3/Scalar-Transport-Equation-Solver/blob/master/documentation/source.JPG "Reaction surce term")

The **initial condition** for the PDE is **_φ(t=0)=0_**.

The **_convective fluxes_** using a third order **QUICK scheme** are given by:

![alt_text](https://github.com/aishkrish3/Scalar-Transport-Equation-Solver/blob/master/documentation/QUICK.JPG "QUICK Scheme")

QUICK (Quadratic Upwind Interpolation for Convective Kinematics) is a simple scheme that is more accurate than a 2nd order finite volume method. It uses piecewise quadratic reconstruction and provides 3rd order accuracy for reconstruction even on non-uniform meshes.

The **_diffusive operator_** is given by:

![alt_text](https://github.com/aishkrish3/Scalar-Transport-Equation-Solver/blob/master/documentation/diffusion.JPG "Diffusive operator")

A 2nd order estimate of derivative over a unifrorm mesh is implemented using a **linear piecewise reconstruction**.

Next, for the **_temporal term_**, denoting the right-hand side of the below equations at timestep n as RHS<sup>n</sup>, the temporal discretization for each scheme is given by:

**First order explicit**:

![alt_text](https://github.com/aishkrish3/Scalar-Transport-Equation-Solver/blob/master/documentation/FirstOrderExplicit.JPG "First Order Explicit")

The first order explicit method implemented is the FTCS (Forward Time Central Space) which is 1st order accurate in time and 2nd order accurate in space. Its stability is kept in check by the Courant–Friedrichs–Lewy condition.

**Crank-Nicolson**:

![alt_text](https://github.com/aishkrish3/Scalar-Transport-Equation-Solver/blob/master/documentation/CrankNicolson.JPG "Crank Nicolson")

The Crank Nicolson can be thought of as an average between FTCS and BTCS (Backward Time Central Space). It is implemented as a two step update as shown above. The first step is explicit and the second step is implicit. It is 2nd order accurate in time and 2nd order accurate in space. The advantage of using Crank Nicolson is that it is always stable and does not depend on the Courant–Friedrichs–Lewy number.


## License:
This is distributed under [MIT License](https://github.com/aishkrish3/Scalar-Transport-Equation-Solver/edit/master/LICENSE).
 
## Author:
Aishwarya Krishnan

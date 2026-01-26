## sphereHPDE_solver
sphereHPDE_solver is a MATLAB code that numerically solves the heat diffusion 
equation in spherical geometry under homogeneous volumetric heating. The solver 
implements a Crank–Nicolson finite-difference scheme with explicit treatment 
of the origin (ghost node) and enforces continuity of both temperature and heat 
flux across the particle–medium interface.

The present v0.1.0-alpha is a preliminary release. It currently addresses the 
case of a single metallic or dielectric sphere embedded in an infinite medium,
subject to a uniform internal heat source Q. Two types of external boundary 
conditions are implemented at the edge of the computational domain:

    +Dirichlet boundary condition: fixed ambient temperature.
    +Neumann boundary condition: zero heat flux (thermally insulated boundary).

This release should be considered work in progress. Future versions will expand
the solver to include further functionalities. Output includes the full spatio-temporal
temperature field inside the sphere and its surrounding medium.


## Project status

This repository is under active development.
Releases are published periodically and archived in Zenodo.


## Citation

This software has been developed to support ongoing research on heat 
diffusion and thermoplasmonic effects in spherical geometries.

A formal publication describing the underlying physical model and
numerical implementation is currently in preparation. you use this 
code in academic work prior to the publication of the corresponding article,
please acknowledge the software by citing the Zenodo record associated with 
the specific release used.



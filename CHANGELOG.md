## v1.1.0 - 2026-10-04

### Added

- Unified solver interface for continuous-wave (CW) and transient thermal calculations.
- Direct stationary formulation for CW heat-diffusion problems.
- Spectral–thermal factorization for efficient reconstruction of wavelength-dependent temperature solutions.
- `temperatureProfile` utility for reconstructing radial temperature profiles at selected wavelengths and times.
- Sparse truncated-inverse preconditioner for improving the conditioning of the thermal system matrices.
- Automated verification suite for testing numerical and physical consistency.
    + Verification of transient-to-CW convergence under constant illumination.
    + Verification against analytical CW thermal solutions.
    + Verification of the quasi-static optical limit.
    + Verification of thermal linearity.
    + Verification of preconditioner performance and preservation of the numerical solution.
- CW example for a gold nanosphere in water comparing perfect thermal contact and finite interfacial thermal conductance.

### Improved

- Separation between optical and thermal contributions to the temperature response.
- Reuse of thermal solutions across multiple wavelengths through spectral factorization.
- Organization of CW and transient calculations under a common solver workflow.
- Numerical verification and reproducibility of the thermal solver.
- Solution-analysis workflow through dedicated post-processing tools.
- Documentation and examples describing the updated solver architecture.

### Changed

- Thermal solutions are now represented in factorized form using optical and boundary contributions instead of requiring explicit reconstruction of the complete wavelength-dependent temperature field.
- CW and transient simulations are now accessed through the common `heatEqSolver` interface using the solver mode configuration.
- Temperature profiles are reconstructed on demand from the factorized solution structure.
- Preconditioning is now available as an optional configurable solver procedure.

## v1.0.0 - 2026-10-01

### Added

- Mie-theory calculation of the optical absorption cross-section using MatScat.
- Support for wavelength-dependent optical properties loaded from material data files.
- Finite interfacial thermal conductance (ITC) between the nanoparticle and the surrounding medium.
- Quasi-static approximation for comparison with the full Mie-theory absorption calculation.
- Validation example comparing the numerical thermal solution with the analytical quasi-static prediction for the maximum temperature increase.
- Repository initialization routine to automatically configure the required MATLAB paths.
- Example script illustrating the definition of the simulation parameters and execution of the solver.

### Improved

- Matrix conditioning procedure used by the Crank–Nicolson thermal solver.
- Organization of the source code into dedicated thermal, optical, tools, examples, and external-dependency directories.
- Handling and interpolation of wavelength-dependent optical material data.
- Code structure and documentation to improve usability and reproducibility.

### Changed

- Optical heating is now determined from the Mie absorption cross-section instead of requiring the absorbed power to be prescribed directly.
- MatScat is included as a third-party dependency for Mie-theory calculations.
- Simulation outputs and examples have been reorganized to provide a clearer workflow for homogeneous spherical nanoparticles.

## v0.1.0 - 2026-01-26

- This release evolved from an earlier preliminary version archived on Zenodo on 2025-09-30.
- Preliminary-version DOI: https://doi.org/10.5281/zenodo.17236155



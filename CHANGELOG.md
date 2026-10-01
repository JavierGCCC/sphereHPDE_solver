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

## v0.1.0 - 2026-01-26 (Original version: Zenodo 2025-09-30)

- Original version. Uploaded to Zenodo on 2025-09-30. DOI: [https://doi.org/10.5281/zenodo.17236155](https://doi.org/10.5281/zenodo.17236155)



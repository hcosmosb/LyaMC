# LyaMC

This code performs a Monte Carlo simulation of Lyα photons in 3D space. The code is entirely written in Fortran, with OpenMP parallelism applied.

This repository provides a vanilla version, which runs with a simplified setup: uniform H I density, temperature, and velocity. Interested users are encouraged to modify the code to read in custom HI density, temperature, and velocity fields.

The detailed description of the code can be found in Park et al. 2022 (https://ui.adsabs.harvard.edu/abs/2022ApJ...931..126P/abstract).

## Repository Contents

[1] run.sh: compilation and execution script

- If using a Fortran compiler other than gfortran, edit this file accordingly.

- Set the OMP_NUM_THREADS variable properly to utilize OpenMP parallelism.

[2] main.f90: main driver code

- Setting PO_plane_flag to True makes the code to save a Lya intensity mapping using the Peeling-off technique.

[3] initial_setup.f90: initial setup for 3D HI density/temperature/velocity grid

- Currently, set to generate uniform HI density of 1e-8 cm^-3, uniform temperature of 2e4 K, and zero velocity.

- Edit this file properly to use custom grid. 

[4] run.f90: main simulation run

- Generates initial photon positions and directions

- Simulates multiple photons (currently, 1e6).

[5] LyART.f90: Monte Carlo Radiative Transfer algorithm

- Core of the simulation

- Several booleans (print_flag, print_flag2, print_flag3) can be set to True for debugging outputs.

[6] rng.f90: random number generator

[7] PO_plane.f90: Peeling off algorithm

[8] exit_condition.f90: condition for finishing run for each photon

[9] Hyojeong.f90, Hyojeong_mod.f90, & CalcT_sub.f90: physical quantity calculation modules written and tested by Hyojeong Kim.



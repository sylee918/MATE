# MATE
Model for Atmospheric Terrestrial Exosphere (MATE) simulation codes are shared here.

# Basic environment
 - MATE is developed in NAS Pleiades system. Please load the following modules.
   - 1) comp-intel/2020.4.304   2) mpi-hpe/mpt.2.28_25Apr23_rhel87
     - Intel compiler is recommended, or you can revise the Makefile for other compilers.
     - MATE is parallelized in MPI. OpenMP is not used.
 - Required input files
    1) indices_for_MSIS_1963-2023_2.txt: Daily OMNI data for Lyman-alpha and photoionization rates (updated by SYLee)
    2) exobase boundary conditions (BC): Daily input files in bindary format are required in current version.
 - Current version does not support running across year. Please run within a time range in a calendar year.


# Update Notes
(09.05.2024) SYLee created this git repository.
 Version: Earth gravity + solar radiation pressure + Coriolis force around the Sun + Photoionization







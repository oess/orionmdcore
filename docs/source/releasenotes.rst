#############
Release Notes
#############


v 1.1.6 June 2022
=================

General Notes
---------------------------------------------------------------------------------
* Add new cube for checking the wall clock time of production runs
* G4/G5 AWS instances support for OpenMM/Gromacs
* Minor bug fixing

New Cubes
--------------------------------------------------------------------------------
* MDSpeedChecker cube. This Cube estimates the wall clock time a production run will 
  require and if the estimated time exceeds the 12 hour time limit, the cube will fail the job. 

Cube Updates
--------------------------------------------------------------------------------
* All cubes have been updated for improved exception handling.

======================


v 1.1.5 April 2022
=================

General Notes
---------------------------------------------------------------------------------
* New cubes for Long Trajectory MD support; these are cubes to enable iterative floes.
* Package dependencies update
* Minor bug fixing

New Cubes
--------------------------------------------------------------------------------
* MDProxyCube. The MD Proxy cube runs chunks of MD simulations in cycles till the desired total MD running time has been reached
* MDRecoveryRestartProxyCube. This cube is able to restart or recover a previously run long MD simulation. It is used in association with the MD proxy cube.

Cube Updates
--------------------------------------------------------------------------------
* All the MD cubes are by default pinning g3 spot instances to improve reliability.

======================


v 1.1.2 December 2021
=================

General Notice
---------------------------------------------------------------------------------
* Better Implicit Solvent Support
* Add CPU support to OpenMM and Gromax 
* Add HMR support to Gromax 

New Cubes
--------------------------------------------------------------------------------
* MDAPIDatasetConverterCube

Cube Updates
--------------------------------------------------------------------------------
* All the MD cubes have been updated to support Implicit Solvent simulations

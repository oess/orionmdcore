#############
Release Notes
#############


v 2.0.0 April 2023
=================

General Notes
---------------------------------------------------------------------------------
* 

New Cubes
--------------------------------------------------------------------------------
* Iterative NPT Cube. This Cube performs NPT MD simulations in the cyclic MD Floe (Plain MD Floe)
  on the provided system. If the MD cycle is complete, the Cube emits the record to the "success" port and 
  if not, it emits the record to the "cycle_out_port" port. 
* Iteration Checker Cube. This Cube checks if the MD cycle is complete. If the MD cycle is complete, 
  it emits the record to the "true" output port and if not, it emits the record 
  to the "false" output port. 

Cube Updates
--------------------------------------------------------------------------------
*  

======================


v 1.1.6 June 2022
=================

General Notes
---------------------------------------------------------------------------------
* Add new Cube for checking the wall clock time of production runs
* G4/G5 AWS instances support for OpenMM/GROMACS
* Minor bug fixing

New Cubes
--------------------------------------------------------------------------------
* MDSpeedChecker Cube. This Cube estimates the wall clock time a production run will 
  require and if the estimated time exceeds the 12 hour time limit, the Cube will fail the job. 

Cube Updates
--------------------------------------------------------------------------------
* All Cubes have been updated for improved exception handling.

======================


v 1.1.5 April 2022
=================

General Notes
---------------------------------------------------------------------------------
* New Cubes for Long Trajectory MD support; these are Cubes to enable iterative Floes.
* Package dependencies update
* Minor bug fixing

New Cubes
--------------------------------------------------------------------------------
* MDProxyCube. The MD Proxy Cube runs chunks of MD simulations in cycles till the desired total MD running time has been reached
* MDRecoveryRestartProxyCube. This Cube is able to restart or recover a previously run long MD simulation. It is used in association with the MD proxy cube.

Cube Updates
--------------------------------------------------------------------------------
* All the MD Cubes are by default pinning g3 spot instances to improve reliability.

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

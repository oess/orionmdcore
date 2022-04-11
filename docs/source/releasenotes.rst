#############
Release Notes
#############


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

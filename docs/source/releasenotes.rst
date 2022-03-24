#############
Release Notes
#############


v 1.1.5 March 2022
=================

General Notice
---------------------------------------------------------------------------------
* New cubes for Long Trajectory MD supports
* Package dependencies update
* Minor bugs fixing

New Cubes
--------------------------------------------------------------------------------
* MDProxyCube. The MD Proxy cube runs chucks of MD simulations in cycles till the desired total md running time has been reached
* MDRecoveryRestartProxyCube. This cube is able to restart or recover a previously run long MD simulation and it is used as support of the MD proxy cube

Cube Updates
--------------------------------------------------------------------------------
* All the MD cubes are by default pinning g3 spot instances

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

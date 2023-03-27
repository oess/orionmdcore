#############
Release Notes
#############

v 2.0.1 March 2023
==================

General Notes
---------------------------------------------------------------------------------
* Since the first MD pkg release we have been using Pickle to serialize and deserialize MD custom objects
  passed along the OpenEye data records. A main side effects of pickling them is that
  the class objects API is embedded inside the produced datasets, therefore updating the object APIs,
  could generate issues. To overcome the problem, all the MD custom objects are now using
  JSON serialization and deserialization which does not have this side effect. To avoid
  compatibility problems, in this release  both serializers (pickle and JSON) are supported however,
  the support to pickle will be dropped in future releases and only the JSON serializer will be kept.
  The MDAPIDatasetConverterCube can be used to convert old datasets into new datasets where the
  JSON serializer/deserializer is used.

New Cubes
--------------------------------------------------------------------------------


Cube Updates
--------------------------------------------------------------------------------
* The MDAPIDatasetConverterCube has been updated to support the conversion from old
  pickle datasets to the new JSON serializer/deserializer
* The ComplexPrepCube has been updated to raise a fatal error if the protein cannot
  been found.

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

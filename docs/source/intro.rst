#############
Introduction
#############

The main idea behind this package to make make it easier for Orion floe developers
to create their own specialized MD floes by providing the core building blocks for
setting up and running MD in Orion, and the data structure for streaming
the results between cubes.
The developers can then focus on their own specialized
application within their own Orion floe package,
which simply installs **orionmdcore** as part of its environment.

The **orionmdcore** package is made available as a conda-installable package
from channel OpenEye/label/Orion; python installation is easy with:
::
    conda install -c conda-forge -c OpenEye/label/Orion -c omnia orionmdcore

However, to run floes built from **orionmdcore** cubes additionally requires the manual
installation of the OpenEye **orionplatform** and **snowball** packages in the conda environment.
Customers that have access to OpenEye MagPie repository and
correctly set their local pip credentials can install them with:
::
     pip install OpenEye-orionplatform[artemis]==4.3.0 OpenEye-toolkits==2021.2.0 OpenEye-snowball==0.23.1

In case of issues please contact OpenEye support.




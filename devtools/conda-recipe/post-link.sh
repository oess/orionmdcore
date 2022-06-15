#!/bin/bash

cat << EOF >> ${PREFIX}/.messages.txt


**************************************************************************
WARNING: The ORIONMDCORE package requires the manual installation of the
OpenEye orionplatform and snowball packages in their working conda environment.

Customers that have access to OpenEye MagPie repository and correctly set
their local pip credentials can have them installed by typing:

pip install OpenEye-orionplatform[artemis]==4.4.0 OpenEye-toolkits==2022.1.0 OpenEye-snowball==0.24.0

In case of issues please contact the OpenEye support.
***************************************************************************
EOF

#!/bin/bash

cat << EOF >> ${PREFIX}/.messages.txt


**************************************************************************
WARNING: The ORIOMDCORE package requires to work the manual installation
of the OpenEye orionplatform and snowball packages in their working conda
environment.

Customers that have access to OpenEye MagPie repository and correctly set
their local pip credentials can have them installed by typing:

pip install OpenEye-orionplatform[artemis]==4.3.0 OpenEye-toolkits==2021.2.0 OpenEye-snowball==0.23.1

In case of issues please contact the OpenEye support
***************************************************************************
EOF

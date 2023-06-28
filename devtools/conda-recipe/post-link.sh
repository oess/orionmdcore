#!/bin/bash

cat << EOF >> ${PREFIX}/.messages.txt


**************************************************************************
WARNING: The ORIONMDCORE package requires the manual installation of the
OpenEye orionplatform and snowball packages in their working conda environment.

Customers that have access to OpenEye MagPie repository and correctly set
their local pip credentials can have them installed by typing:

pip install "openeye-orionplatform[artemis]>=5.1.0,<6.0.0" openeye-toolkits>=2023.1.0 Openeye-snowball==0.26.0

In case of issues please contact the OpenEye support.
***************************************************************************
EOF

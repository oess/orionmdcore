#!/bin/bash

cat << EOF >> ${PREFIX}/.messages.txt


**************************************************************************
WARNING: The ORIOMDCORE package requires to work the manual installation
of the OpenEye orionplatform and snowball packages.

Customers that have access and correctly set their local pip credentials
can have them installed by typing:

pip install OpenEye-orionplatform[artemis]==4.0.0 OpenEye-snowball==0.21.0

In case of issues please contact the OpenEye support
***************************************************************************
EOF
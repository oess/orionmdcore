# (C) 2022 Cadence Design Systems, Inc. (Cadence) 
# All rights reserved.
#
# TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# provided to current licensees or subscribers of Cadence products or
# SaaS offerings (each a "Customer").
# Customer is hereby permitted to use, copy, and modify the Sample Code,
# subject to these terms. Cadence claims no rights to Customer's
# modifications. Modification of Sample Code is at Customer's sole and
# exclusive risk. Sample Code may require Customer to have a then
# current license or subscription to the applicable Cadence offering.
# THE SAMPLE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED.  OPENEYE DISCLAIMS ALL WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. In no event shall Cadence be
# liable for any damages or liability in connection with the Sample Code
# or its use.

__author__ = "OpenEye, Cadence Molecular Sciences"
__email__ = "support@eyesopen.com"
__version__ = '2.0.0a7'

__installation__error__ = """
The ORIONMDCORE package requires the manual installation of the
OpenEye orionplatform and snowball packages in their working conda environment.

Customers that have access to OpenEye MagPie repository and correctly set
their local pip credentials can have them installed by typing:

pip install OpenEye-orionplatform[artemis]==4.4.0 OpenEye-toolkits==2022.1.1 OpenEye-snowball==0.24.0

In case of issues please contact the OpenEye support.
"""

# try:
#     import openeye
#     import snowball
#     import orionplatform
# except ImportError:
#     raise ImportError(__installation__error__)



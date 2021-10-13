# (C) 2021 OpenEye Scientific Software Inc. All rights reserved.
#
# TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# provided to current licensees or subscribers of OpenEye products or
# SaaS offerings (each a "Customer").
# Customer is hereby permitted to use, copy, and modify the Sample Code,
# subject to these terms. OpenEye claims no rights to Customer's
# modifications. Modification of Sample Code is at Customer's sole and
# exclusive risk. Sample Code may require Customer to have a then
# current license or subscription to the applicable OpenEye offering.
# THE SAMPLE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED.  OPENEYE DISCLAIMS ALL WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. In no event shall OpenEye be
# liable for any damages or liability in connection with the Sample Code
# or its use.


__author__ = "Gaetano Calabro"
__email__ = "gcalabro@eyesopen.com"
__version__ = '1.1.0b6'

__installation__error__ = """
ERROR: The ORIOMDCORE package requires to work the manual installation of the OpenEye 
orionplatform and snowball packages. Customers that have access and correctly set
their local pip credentials can have them installed typing:

pip install OpenEye-orionplatform[artemis]==4.0.0 OpenEye-snowball==0.21.0

In case of issues please contact the OpenEye support
"""

# try:
#     import openeye
#     import snowball
#     import orionplatform
# except ImportError:
#     raise ImportError(__installation__error__)

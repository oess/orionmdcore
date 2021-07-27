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

import sys

from orionmdcore import __version__

from setuptools import setup, find_packages

if sys.argv[-1] == 'setup.py':
    print("To install, run 'python setup.py install'")
    print()

if sys.version_info[:3] < (3, 0):
    print("OEOMMTools requires Python 3.0 or later (%d.%d detected)." %
          sys.version_info[:2])
    sys.exit(-1)


descr = """
Core functionalities for MD developing in OpenEye Orion 
"""

setup(
    name                 ='orionmdcore',
    version              =__version__,
    description          ='OpenEye Orion Molecular Dynamics Core',
    long_description     =descr,
    url                  ='https://github.com/oess/orionmdcore',
    author               ='Gaetano Calabro and Christopher Bayly',
    author_email         ='gcalabro -at- eyesopen.com',
    platforms            =['Linux-64', 'Mac OSX-64', 'Unix-64'],
    packages             =find_packages()+['tests'],
    include_package_data = True,
    zip_safe             = False
)

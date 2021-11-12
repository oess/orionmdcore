# (C) 2020 OpenEye Scientific Software Inc. All rights reserved.
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


import simtk

from simtk import unit

import itertools

md_keys_converter = {'OpenMM':

                         {'constraints':
                              {'None': 'None', 'Bonds2H': 'HBonds', 'Angles2H': 'HAngles', 'All-Bonds': 'AllBonds'}
                          },

                     'Gromacs':
                         {'constraints':
                              {'None': 'none', 'Bonds2H': 'h-bonds', 'Angles2H': 'h-angles', 'All-Bonds': 'all-bonds'}
                          }
                     }


class MDState(object):

    def __init__(self, parmed_structure):

        if not parmed_structure.positions:
            raise RuntimeError('Atom positions are not defined')
        else:
            # The returned object is an OpenMM Quantity with units
            self.__positions__ = parmed_structure.positions

        if parmed_structure.velocities is None:
            self.__velocities__ = None
        else:
            # Parmed stores the velocities as a numpy array in unit of angstrom/picoseconds
            self.__velocities__ = parmed_structure.velocities * simtk.unit.angstrom/simtk.unit.picosecond
            # The returned object is an OpenMM Quantity with units

        if parmed_structure.box_vectors is None:
            self.__box_vectors__ = None
        else:
            self.__box_vectors__ = parmed_structure.box_vectors

    def get_positions(self):
        return self.__positions__

    def get_oe_positions(self):
        pos = self.__positions__.in_units_of(unit.angstrom) / unit.angstrom
        pos = list(itertools.chain.from_iterable(pos))
        return pos

    def get_velocities(self):
        return self.__velocities__

    def get_box_vectors(self):
        return self.__box_vectors__

    def set_positions(self, positions):
        if isinstance(positions, simtk.unit.quantity.Quantity):
            self.__positions__ = positions
        else:
            raise ValueError("It was not possible to set the positions")

        return

    def set_velocities(self, velocities):
        if isinstance(velocities, simtk.unit.quantity.Quantity):
            self.__velocities__ = velocities
        else:
            raise ValueError("It was not possible to set the velocities")

        return

    def set_box_vectors(self, box_vectors):
        if isinstance(box_vectors, simtk.unit.quantity.Quantity):
            self.__box_vectors__ = box_vectors
        else:
            raise ValueError("It was not possible to set the box vectors")

        return

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
from abc import ABC, abstractmethod

import parmed

from simtk import unit

from simtk.openmm import Vec3

import itertools

from orionmdcore.mdengine import json_numpy

import numpy as np


class MDState(object):
    def __init__(self, parmed_structure=None):

        if parmed_structure is not None:

            if parmed_structure.coordinates is None:
                raise RuntimeError("Atom positions are not defined")
            else:  # numpy array in Angstroms unit
                self.__coordinates__ = parmed_structure.coordinates

            if parmed_structure.velocities is None:
                self.__velocities__ = None
            else:
                # Parmed stores the velocities as a numpy array in unit of angstrom/picoseconds
                self.__velocities__ = parmed_structure.velocities

            if parmed_structure.box is None:
                self.__box__ = None
            else:
                # Parmed Numpy Array [lx, ly, lz, 90, 90, 90]
                self.__box__ = parmed_structure.box
        else:
            return

    def __getstate__(self):
        state = dict()

        state["coordinates"] = json_numpy.dumps(self.__coordinates__)

        state["box"] = (
            json_numpy.dumps(self.__box__) if self.__box__ is not None else None
        )

        state["velocities"] = (
            json_numpy.dumps(self.__velocities__) if self.__velocities__ is not None else None
        )

        return state

    def __copy__(self):
        new_md_state = MDState()
        new_md_state.__coordinates__ = self.__coordinates__
        new_md_state.__velocities__ = self.__velocities__
        new_md_state.__box__ = self.__box__

        return new_md_state

    def __setstate_new_api__(self, state):
        for state_name, state_val in state.items():

            if state_name == "coordinates":
                coordinates = json_numpy.loads(state_val)
                self.__coordinates__ = coordinates

            if state_name == "box":
                if state_val is not None:
                    box = json_numpy.loads(state_val)
                    self.__box__ = box
                else:
                    self.__box__ = None

            if state_name == "velocities":
                if state_val is not None:
                    velocities = json_numpy.loads(state_val)
                    self.__velocities__ = velocities
                else:
                    self.__velocities__ = None

    def get_positions(self):

        # OLD API TRAP
        if '__positions__' in self.__dict__:
            return self.__positions__

        elif self.__coordinates__ is not None:
            return [Vec3(coord[0], coord[1], coord[2]) for coord in self.__coordinates__] * unit.angstroms
        else:
            return None

    def get_oe_positions(self):

        # OLD API TRAP
        if '__positions__' in self.__dict__:
            return list(itertools.chain.from_iterable(self.__positions__.value_in_unit(unit.angstrom)))

        if self.__coordinates__ is not None:
            return list(itertools.chain.from_iterable(self.__coordinates__))
        else:
            return None

    def get_velocities(self):

        # OLD API TRAP
        if unit.is_quantity(self.__velocities__):
            return self.__velocities__

        if self.__velocities__ is not None:
            return self.__velocities__ * unit.angstrom / unit.picosecond
        else:
            return None

    def get_box_vectors(self):

        # OLD API TRAP
        if '__box_vectors__' in self.__dict__:
            return self.__box_vectors__

        elif self.__box__ is not None:

            av = [self.__box__[0], 0.0, 0.0]
            bv = [0.0, self.__box__[1], 0.0]
            cv = [0.0, 0.0, self.__box__[2]]

            return (av, bv, cv) * unit.angstrom

        else:
            return None

    def set_positions(self, positions):
        if isinstance(positions, unit.quantity.Quantity):
            self.__coordinates__ = np.array([[coord[0], coord[1], coord[2]] for coord in
                                             positions.value_in_unit(unit.angstrom)])
        else:
            raise ValueError("It was not possible to set the positions")

        return

    def set_velocities(self, velocities):
        if isinstance(velocities, unit.quantity.Quantity):
            self.__velocities__ = np.array([[vel[0], vel[1], vel[2]] for vel in
                                            velocities.value_in_unit(unit.angstrom / unit.picoseconds)])
        else:
            raise ValueError("It was not possible to set the velocities")

        return

    def set_box_vectors(self, box_vectors):
        if isinstance(box_vectors, unit.quantity.Quantity):
            self.__box__ = np.array([box_vectors[0][0].value_in_unit(unit.angstrom),
                                     box_vectors[1][1].value_in_unit(unit.angstrom),
                                     box_vectors[2][2].value_in_unit(unit.angstrom), 90.0, 90.0, 90.0],
                                    dtype=np.float64)
        else:
            raise ValueError("It was not possible to set the box vectors")

        return


class MDSimulations(ABC):
    @abstractmethod
    def __init__(self, mdstate, ff_parameters, opt):

        if not isinstance(mdstate, MDState):

            if str(type(mdstate)) == '<class \'MDOrion.MDEngines.utils.MDState\'>':

                print("\n\nDEPRECATION WARNING:\n"
                      "You are trying to use a dataset produced by using the old openmm_orion API "
                      "which will not be supported in future releases. "
                      "The API has been moved to the new orionmdcore pkg: "
                      "https://github.com/oess/orionmdcore\n\n", flush=True)
            else:
                raise ValueError("{} is not a valid MDState Object".format(type(mdstate)))

        if not isinstance(ff_parameters, parmed.Structure):
            raise ValueError(
                "{} is not a Parmed Structure Object".format(type(ff_parameters))
            )

        if not isinstance(opt, dict):
            raise ValueError("{} is not a dictionary".format(type(opt)))

    @abstractmethod
    def run(self):
        pass

    @abstractmethod
    def update_state(self):
        pass

    @abstractmethod
    def clean_up(self):
        pass

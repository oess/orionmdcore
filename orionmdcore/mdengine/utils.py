from abc import ABC, abstractmethod

import parmed

from simtk.openmm import Vec3

from simtk import unit

import itertools

import numpy as np

# import json_numpy
from orionmdcore.mdengine import json_numpy


class MDState(object):

    def __init__(self, from_parmed=None):

        if from_parmed is not None:
            if from_parmed.coordinates is None:
                raise RuntimeError("Atom positions are not defined")
            else:
                # Coordinate in Angstrom
                self.__coordinates__ = from_parmed.coordinates

            if from_parmed.velocities is None:
                self.__velocities__ = None
            else:
                # Parmed stores the velocities as a numpy array in unit of angstrom/picoseconds
                self.__velocities__ = from_parmed.velocities

            if from_parmed.box is None:
                # Lengths in Angstroms
                self.__box__ = None
            else:
                self.__box__ = np.array([from_parmed.box[0],
                                         from_parmed.box[1],
                                         from_parmed.box[2]],
                                        dtype=np.float64)

        else:
            return

    def __getstate__(self):
        state = dict()

        state["positions"] = json_numpy.dumps(self.__coordinates__)

        if self.__box__ is not None:
            state["box_vectors"] = json_numpy.dumps(self.__box__)
        else:
            state["box_vectors"] = None

        if self.__velocities__ is not None:

            state["velocities"] = json_numpy.dumps(self.__velocities__)
        else:
            state["velocities"] = None

        return state

    def __setstate__(self, state):
        for state_name, state_val in state.items():

            if state_name == "positions":
                positions = json_numpy.loads(state_val)
                self.__coordinates__ = positions

            if state_name == "box_vectors":
                if state_val is not None:
                    box_vec = json_numpy.loads(state_val)
                    self.__box__ = box_vec
                else:
                    self.__box_vectors__ = None

            if state_name == "velocities":
                if state_val is not None:
                    velocities = json_numpy.loads(state_val)
                    self.__velocities__ = velocities
                else:
                    self.__velocities__ = None

    def get_positions(self):
        if self.__coordinates__ is not None:
            return [Vec3(coord[0], coord[1], coord[2]) for coord in self.__coordinates__] * unit.angstroms
        else:
            return None

    def get_oe_positions(self):
        if self.__coordinates__ is not None:
            return list(itertools.chain.from_iterable(self.__coordinates__))
        else:
            return None

    def get_velocities(self):
        if self.__velocities__ is not None:
            return self.__velocities__ * unit.angstrom / unit.picosecond
        else:
            return None

    def get_box_vectors(self):

        if self.__box__ is not None:

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
                                     box_vectors[2][2].value_in_unit(unit.angstrom)],
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

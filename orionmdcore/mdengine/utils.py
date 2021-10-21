from abc import ABC, abstractmethod

import parmed

import simtk

from simtk import unit

import itertools


class MDState(object):
    def __init__(self, parmed_structure):

        if not parmed_structure.positions:
            raise RuntimeError("Atom positions are not defined")
        else:
            # The returned object is an OpenMM Quantity with units
            self.__positions__ = parmed_structure.positions

        if parmed_structure.velocities is None:
            self.__velocities__ = None
        else:
            # Parmed stores the velocities as a numpy array in unit of angstrom/picoseconds
            self.__velocities__ = (
                parmed_structure.velocities
                * simtk.unit.angstrom
                / simtk.unit.picosecond
            )
            # The returned object is an OpenMM Quantity with units

        if parmed_structure.box_vectors is None:
            self.__box_vectors__ = None
        else:
            self.__box_vectors__ = parmed_structure.box_vectors

    # def __getstate__(self):
    #     state = dict()
    #
    #     state["positions"] = data_utils.encodePyObj(self.__positions__)
    #
    #     state["box_vectors"] = (
    #         data_utils.encodePyObj(self.__box_vectors__)
    #         if self.__box_vectors__
    #         else None
    #     )
    #
    #     state["velocities"] = (
    #         data_utils.encodePyObj(self.__velocities__) if self.__velocities__ else None
    #     )
    #
    #     return state
    #
    # def __setstate__(self, state):
    #     for state_name, state_val in state.items():
    #
    #         if state_name == "positions":
    #             positions = data_utils.decodePyObj(state_val)
    #             self.__positions__ = positions
    #
    #         if state_name == "box_vectors":
    #             if state_val is not None:
    #                 box_vec = data_utils.decodePyObj(state_val)
    #                 self.__box_vectors__ = box_vec
    #             else:
    #                 self.__box_vectors__ = None
    #
    #         if state_name == "velocities":
    #             if state_val is not None:
    #                 velocities = data_utils.decodePyObj(state_val)
    #                 self.__velocities__ = velocities
    #             else:
    #                 self.__velocities__ = None

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


class MDSimulations(ABC):
    @abstractmethod
    def __init__(self, mdstate, ff_parameters, opt):

        if not isinstance(mdstate, MDState):
            raise ValueError("{} is not a MDState Object".format(type(mdstate)))

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

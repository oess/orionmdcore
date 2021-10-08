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

try:
    from orionclient.session import in_orion
except ImportError:
    from orionmdcore import __installation__error__
    raise ImportError(__installation__error__)

import time

import fcntl

import os

from simtk import unit

import math


md_keys_converter = {
    "OpenMM": {
        "constraints": {
            "None": "None",
            "Bonds2H": "HBonds",
            "Angles2H": "HAngles",
            "All-Bonds": "AllBonds",
        }
    },
    "Gromacs": {
        "constraints": {
            "None": "none",
            "Bonds2H": "h-bonds",
            "Angles2H": "h-angles",
            "All-Bonds": "all-bonds",
        }
    },
}


def local_cluster(sim):
    def wrapper(*args):

        mdstate = args[0]
        ff_parameters = args[1]
        opt = args[2]

        if "OE_VISIBLE_DEVICES" in os.environ and not in_orion():

            gpus_available_indexes = os.environ["OE_VISIBLE_DEVICES"].split(",")

            opt["Logger"].info("OE LOCAL FLOE CLUSTER OPTION IN USE")

            if "OE_MAX" in os.environ:
                opt["OE_MAX"] = int(os.environ["OE_MAX"])
            else:
                opt["OE_MAX"] = 1

            opt["Logger"].info("OE MAX = {}".format(opt["OE_MAX"]))

            while True:

                for gpu_id in gpus_available_indexes:

                    for p in range(0, opt["OE_MAX"]):

                        fn = str(gpu_id) + "_" + str(p) + ".txt"

                        try:
                            with open(fn, "a") as file:

                                fcntl.flock(file, fcntl.LOCK_EX | fcntl.LOCK_NB)
                                # opt['Logger'].warn("LOCKED GPU ID = {} - MOL ID = {}".format(gpu_id, opt['system_id']))
                                file.write(
                                    "MD - name = {} MOL_ID = {} GPU_IDS = {} GPU_ID = {}\n".format(
                                        opt["system_title"],
                                        opt["system_id"],
                                        gpus_available_indexes,
                                        str(gpu_id),
                                    )
                                )
                                opt["gpu_id"] = str(gpu_id)

                                new_mdstate = sim(mdstate, ff_parameters, opt)

                                time.sleep(5.0)
                                # opt['Logger'].warn("UNLOCKING GPU ID = {} - MOL ID = {}".format(gpu_id, opt['system_id']))
                                fcntl.flock(file, fcntl.LOCK_UN)
                                return new_mdstate

                        except BlockingIOError:
                            time.sleep(0.1)

                        except Exception as e:  # If the simulation fails for other reasons
                            try:
                                time.sleep(5.0)
                                fcntl.flock(file, fcntl.LOCK_UN)
                            except Exception as e:
                                pass
                            raise ValueError("{} Simulation Failed".format(str(e)))
        else:
            new_mdstate = sim(*args)
            return new_mdstate

    return wrapper


@local_cluster
def md_simulation(mdstate, ff_parameters, opt):

    if opt["md_engine"] == "OpenMM":

        from orionmdcore.cubes.md.openmm.simtools import OpenMMSimulations

        MDSim = OpenMMSimulations(mdstate, ff_parameters, opt)

        MDSim.run()

        new_mdstate = MDSim.update_state()

        MDSim.clean_up()

        return new_mdstate

    if opt["md_engine"] == "Gromacs":

        from orionmdcore.cubes.md.gromacs.simtools import GromacsSimulations

        MDSim = GromacsSimulations(mdstate, ff_parameters, opt)

        MDSim.run()

        new_mdstate = MDSim.update_state()

        MDSim.clean_up()

        return new_mdstate

    else:
        raise ValueError(
            "The selected MD engine is not currently supported: {}".format(
                opt["md_engine"]
            )
        )


def update_cube_parameters_in_place(record, parameter_dic):

    from orionmdcore.standards.standards import Fields

    if record.has_value(Fields.cube_parameters_update):

        cube_parameters_dic = record.get_value(Fields.cube_parameters_update)

        parameters_intersection_set = cube_parameters_dic.keys() & parameter_dic.keys()

        for p in parameters_intersection_set:
            parameter_dic["Logger"].info(
                "Updating in place cube parameter: {} from {} to {}".format(
                    p, parameter_dic[p], cube_parameters_dic[p]
                )
            )
        parameter_dic.update(
            {k: cube_parameters_dic[k] for k in parameters_intersection_set}
        )


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
                "Updating cube parameter: {} from {} to {}\n".format(
                    p, parameter_dic[p], cube_parameters_dic[p]
                )
            )
        parameter_dic.update(
            {k: cube_parameters_dic[k] for k in parameters_intersection_set}
        )


def schedule_cycles(cube_param_dic, record_info_dic):

    speed_ns_per_day = record_info_dic["speed_ns_per_day"]
    hmr = record_info_dic["hmr"]
    md_engine = record_info_dic["md_engine"]

    if md_engine == "Gromacs":
        time_step = 2 * unit.femtoseconds
    else:
        if hmr:
            time_step = 4 * unit.femtoseconds
        else:
            time_step = 2 * unit.femtoseconds

    cube_max_run_time = cube_param_dic["cube_max_run_time"] * unit.hours
    time = cube_param_dic["time"] * unit.nanoseconds
    time_ns = time.in_units_of(unit.nanoseconds)
    total_steps = math.ceil(time_ns / time_step)

    trajectory_interval_ns = cube_param_dic["trajectory_interval"] * unit.nanoseconds

    # Per cycle md running time
    running_time_per_cycle = (
        cube_max_run_time.value_in_unit(unit.day) * speed_ns_per_day * unit.nanoseconds
    )
    running_time_per_cycle_in_ns = running_time_per_cycle.in_units_of(unit.nanoseconds)

    md_steps_per_cycle = int(running_time_per_cycle_in_ns / time_step)

    cycles = int(time_ns / running_time_per_cycle_in_ns)
    cycle_rem = time_ns.value_in_unit(
        unit.nanoseconds
    ) % running_time_per_cycle_in_ns.value_in_unit(unit.nanoseconds)

    trajectory_steps = int(
        round(
            trajectory_interval_ns.value_in_unit(unit.nanosecond)
            / (time_step.value_in_unit(unit.nanoseconds))
        )
    )

    schedule = dict()
    if cycles == 0:
        time_step_to_time = (
            int(round(cycle_rem / time_step.value_in_unit(unit.nanoseconds)))
            * time_step
        )

        frame_per_cycle = int(total_steps / trajectory_steps)

        schedule[0] = (
            time_step_to_time.value_in_unit(unit.nanosecond),
            frame_per_cycle,
        )
    else:
        time_steps_to_time = md_steps_per_cycle * time_step.in_units_of(
            unit.nanoseconds
        )
        last_time = (total_steps - md_steps_per_cycle * cycles) * time_step.in_units_of(
            unit.nanoseconds
        )

        frames_per_cycle = int(md_steps_per_cycle / trajectory_steps)
        frame_per_cycle_dic = {i: frames_per_cycle for i in range(0, cycles)}
        frame_per_cycle_dic[cycles] = int(
            (total_steps - md_steps_per_cycle * cycles) / trajectory_steps
        )

        schedule = {
            i: (
                time_steps_to_time.value_in_unit(unit.nanosecond),
                frame_per_cycle_dic[i],
            )
            for i in range(0, cycles)
        }
        schedule[cycles] = (
            last_time.value_in_unit(unit.nanosecond),
            frame_per_cycle_dic[cycles],
        )

    return schedule

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

try:
    from openeye import oechem

    from orionplatform.mixins import RecordPortsMixin

    from floe.api import ParallelMixin, parameters, ComputeCube

    from orionplatform.ports import RecordInputPort, RecordOutputPort

    from snowball.utils.log_params import LogFieldParam

except ImportError:
    from orionmdcore import __installation__error__
    raise ImportError(__installation__error__)

import traceback

from orionmdcore.standards.standards import (
    MDStageTypes,
    MDEngines,
    MDStageNames,
    Fields,
    MDFileNames
)

from orionmdcore.mdrecord.mdrecord import MDDataRecord

from orionmdcore.cubes.md.utils import (
    md_simulation,
    update_cube_parameters_in_place,
)

import copy

import textwrap

import os

from orionmdcore.standards.utils import check_filename

from orionmdcore.cubes.md.utils import str_schedule, schedule_cycles, size_dependent_schedule_cycles

from snowball import GPU_INSTANCE_TYPES_SUPPORT_FOR_MD

from orionclient.session import in_orion

import mdtraj

import glob

import gzip

import tarfile

from tempfile import TemporaryDirectory

import json

import tempfile

import h5py


class IterMDNptCube(RecordPortsMixin, ComputeCube):
    title = "Iterative NPT Cube"

    classification = [["MD Simulations"]]
    tags = ["Gromacs", "OpenMM", "NPT"]

    description = """
    This Cube performs NPT MD simulations in the cyclic MD Floe (Plain MD Floe)
    on the provided system.
    The system must have been parametrized by the Force Field cube and the system Parmed 
    object must be present on the input record. In addition, a system identification 
    number must be present on the input record as well. This can be accomplished 
    by using the “ID Setting Cube”. The NPT MD simulation is performed by the selected 
    MD engine, currently OpenMM and Gromacs only. Restraints and constraints can be 
    used as well. Currently implicit solvent models can be used in OpenMM only. 
    The cube requires a record as input and produces a new record with the time evolved 
    system. The total sampling time can be set by using the “time” cube parameter 
    while the trajectory snapshots can be set by using the “trajectory_interval” cube 
    parameter.
    If the MD cycle is complete, the Cube emits the record to the "success" port and 
    if not, it emits the record to the "cycle_out_port" port. 
    """

    # for Exception Handler
    log_field = LogFieldParam()

    uuid = "8ea50706-2560-4eb4-a367-cb585e636c84"

    # Override defaults for some parameters
    parameter_overrides = {
        "gpu_count": {"default": 1},
        "cpu_count": {"default": 16},
        "instance_type": {"default": GPU_INSTANCE_TYPES_SUPPORT_FOR_MD},  # Gpu Family selection
        # "instance_type": {"default": ""},  # Gpu Family selection
        "disk_space": {"default": 20000},
        "memory_mb": {"default": 14000},
        "spot_policy": {"default": "Prohibited"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1},  # 1 molecule at a time
    }

    temperature = parameters.DecimalParameter(
        "temperature", default=300.0, help_text="Temperature (Kelvin)"
    )

    pressure = parameters.DecimalParameter(
        "pressure", default=1.0, help_text="Pressure (atm)"
    )

    time = parameters.DecimalParameter(
        "time", default=10, help_text="NPT simulation time in nanoseconds"
    )

    restraints = parameters.StringParameter(
        "restraints",
        default="",
        help_text=""""Mask selection to apply harmonic restraints. 
        Possible keywords are: ligand, protein, water, ions, 
        ca_protein, cofactors. The selection can be refined 
        by using logical tokens: not, noh, and, or, diff, around""",
    )

    restraintWt = parameters.DecimalParameter(
        "restraintWt",
        default=2.0,
        help_text="Restraint weight for xyz atom restraints in kcal/(mol A^2)",
    )

    restraint_to_reference = parameters.BooleanParameter(
        "restraint_to_reference",
        title="Restraint to reference",
        default=True,
        help_text="If True the starting reference system coordinates will be used "
        "to restraint the system",
    )

    freeze = parameters.StringParameter(
        "freeze",
        default="",
        help_text="""Mask selection to freeze atoms along the MD simulation.
        Possible keywords are: ligand, protein, water, ions, ca_protein,
        cofactors. The selection can be refined by using logical tokens:
        not, noh, and, or, diff, around. Not currently implemented in Gromacs""",
    )

    nonbondedCutoff = parameters.DecimalParameter(
        "nonbondedCutoff",
        default=10,
        help_text="""The non-bonded cutoff in angstroms.
        This is ignored if non-bonded method is NoCutoff""",
    )

    constraints = parameters.StringParameter(
        "constraints",
        default="Bonds2H",
        choices=["None", "Bonds2H", "Angles2H", "All-Bonds"],
        help_text="""None, Bonds2H, Angles2H, or All-Bonds
        Which type of constraints to add to the system.
        None means no bonds are constrained.
        HBonds means bonds with hydrogen are constrained, etc.""",
    )

    implicit_solvent = parameters.StringParameter(
        "implicit_solvent",
        title="Implicit Solvent",
        default="None",
        choices=["None", "HCT", "OBC1", "OBC2", "GBn", "GBn2"],
        help_text="Implicit Solvent Model. Not Currently implemented in Gromacs",
    )

    trajectory_interval = parameters.DecimalParameter(
        "trajectory_interval",
        title="Trajectory interval",
        default=0.01,
        help_text="""Time interval for trajectory snapshots in ns. 
        If 0 the trajectory file will not be generated""",
    )

    reporter_interval = parameters.DecimalParameter(
        "reporter_interval",
        title="Reporter interval",
        default=0.01,
        help_text="""Time interval for reporting data in ns. 
        If 0 the reporter file will not be generated""",
    )

    trajectory_frames = parameters.IntegerParameter(
        "trajectory_frames",
        title="Trajectory frames",
        default=0,
        help_text="""The total number of trajectory frames. If it is 
            set to zero and the trajectory interval parameter is set 
            to zero no trajectory is generated. If it is different from zero 
            and the trajectory interval parameter is set to zero the produced 
            trajectory will have the selected number of frames. If different
            from zero and the trajectory interval parameter is different from 
            zero the total number of generated frames will be calculated by just 
            using the trajectory interval and the md time step (2fs and 4fs hmr on)""",
    )

    suffix = parameters.StringParameter(
        "suffix", default="npt", help_text="Filename suffix for output simulation files"
    )

    center = parameters.BooleanParameter(
        "center", default=False, help_text="Center the system to the OpenMM unit cell"
    )

    verbose = parameters.BooleanParameter(
        "verbose", default=True, help_text="Increase log file verbosity"
    )

    hmr = parameters.BooleanParameter(
        "hmr",
        title='Hydrogen Mass Repartitioning',
        default=True,
        help_text="Give hydrogens more mass to speed up the MD",
    )

    save_md_stage = parameters.BooleanParameter(
        "save_md_stage",
        title="Save MD stage",
        default=False,
        help_text="""Save the MD simulation stage. If True the MD,
           simulation data will be appended to the md simulation stages 
           otherwise the last MD stage will be overwritten""",
    )

    md_engine = parameters.StringParameter(
        "md_engine",
        title="MD engine",
        default="OpenMM",
        choices=["OpenMM", "Gromacs"],
        help_text="Select the available MD engine",
    )

    cube_max_run_time = parameters.DecimalParameter(
        "cube_max_run_time", title='Cube max run time', default=1, help_text="Max Cube Running Time in hrs. Max 11hrs"
    )

    speed_check_time = parameters.DecimalParameter(
        "speed_check_time", default=0.05, help_text="Speed check simulation time in nanoseconds"
    )

    cycle_out_port = RecordOutputPort(
        "cycle_out_port",
    )

    def begin(self):
        self.opt = vars(self.args)
        self.opt["Logger"] = self.log
        self.opt["SimType"] = "npt"

        return
    
    def process(self, record, port): 
        try:
            opt = dict(self.opt)
            opt["CubeTitle"] = self.title    

            # Update Cube Parameters
            update_cube_parameters_in_place(record, opt)

            # Logger string
            str_logger = "-" * 32 + " NPT CUBE PARAMETERS " + "-" * 32
            str_logger += "\n{:<25} = {:<10}".format("Cube Title", opt["CubeTitle"])

            for k, v in sorted(self.parameters().items()):
                if k != 'log_field':
                    tmp_default = copy.deepcopy(v)

                    if v.default is None:
                        tmp_default.default = "None"
                    elif isinstance(v, parameters.BooleanParameter):
                        if v.default:
                            tmp_default.default = "True"
                        else:
                            tmp_default.default = "False"
                    else:
                        tmp_description = textwrap.fill(
                            " ".join(v.description.split()),
                            subsequent_indent=" " * 39,
                            width=80,
                        )
                        str_logger += "\n{:<25} = {:<10} {}".format(
                            k, getattr(self.args, tmp_default.name), tmp_description
                        )

            str_logger += "\n{:<25} = {:<10}".format("Simulation Type", opt["SimType"])

            # Create the MD record to use the MD Record API
            mdrecord = MDDataRecord(record)

            system_title = mdrecord.get_title

            opt["system_title"] = system_title
            opt["system_id"] = mdrecord.get_flask_id

            flask = mdrecord.get_stage_topology()
            mdstate = mdrecord.get_stage_state()

            if opt["restraint_to_reference"]:
                opt["reference_state"] = mdrecord.get_stage_state(
                    stg_name=MDStageNames.ForceField
                )

            opt["out_directory"] = mdrecord.cwd
            opt["molecule"] = flask
            opt["str_logger"] = str_logger
            opt["Logger"].info(
                "[{}] START TEST NPT SIMULATION: {}".format(opt["CubeTitle"], system_title)
            )

            opt["out_fn"] = check_filename(
                os.path.basename(opt["out_directory"])
                + "_"
                + opt["system_title"]
                + "_"
                + str(opt["system_id"])
                + "-"
                + opt["suffix"]
            )

            # Trajectory file name if any generated
            opt["trj_fn"] = opt["out_fn"] + "_" + "traj.tar.gz"

            # Extract the Parmed structure and synchronize it with the last MD stage state
            parmed_structure = mdrecord.get_parmed(sync_stage_name="last")

            #########################################
            # 1. Run 0.02ns MD to measure the speed #
            #########################################

            # a. Grep instance type
            if in_orion():
                current_instance = os.environ["ORION_INSTANCE"]
            else:
                current_instance = "local"
            # b. Check if the speed of the instance type can be found in instance_speed_dic_fld 
            new_instance = True
            if record.has_value(Fields.instance_speed_dic_fld):
                instance_speed_dic = record.get_value(Fields.instance_speed_dic_fld)
                if current_instance in instance_speed_dic.keys():
                    new_instance = False
            else: 
                instance_speed_dic = dict()
            # c. Check the time remaining 
            time_nonzero = True
            if opt["time"] == 0:
                time_nonzero = False
                if current_instance not in instance_speed_dic.keys():
                    print("speed is unknown, but since no steps left, set speed 1.") ## should change!!
                    instance_speed_dic[current_instance] = 1

            # d. If speed info cant be found & time remaining >0: (1) run 0.02 ns;(2)gen info_dic;(3)save speed 
            if time_nonzero is True and new_instance is True: 
                # (1) run 0.02 ns 
                tmp_opt = dict()
                for k, v in opt.items():
                    tmp_opt[k] = v 

                tmp_opt["time"] = opt["speed_check_time"]
                tmp_opt["suffix"] = "tmp_run"
                new_mdstate = md_simulation(mdstate, parmed_structure, tmp_opt)

                # (2) gen info_dic
                info_dic = dict()
                info_dic["speed_ns_per_day"] = tmp_opt["speed_ns_per_day"]
                info_dic["instance_type"] = tmp_opt["instance_type"]
                info_dic["hmr"] = tmp_opt["hmr"]
                info_dic["md_engine"] = tmp_opt["md_engine"]


                # (3) save speed in instance_speed_dic 
                instance_speed_dic[current_instance] = tmp_opt["speed_ns_per_day"]
                record.set_value(Fields.instance_speed_dic_fld, instance_speed_dic)
            else: 
                info_dic = dict()
                info_dic["speed_ns_per_day"] = instance_speed_dic[current_instance]
                info_dic["instance_type"] = opt["instance_type"]
                info_dic["hmr"] = opt["hmr"]
                info_dic["md_engine"] = opt["md_engine"]
            
            # Save flask size for cube run time
            info_dic["num_atoms"] = flask.NumAtoms()

            #################################
            # 2. Schedule the current cycle #
            #################################

            Terminate=True

            if "speed_ns_per_day" not in info_dic:
                raise ValueError("Last MD stage does not have speed information")

            flask_title = record.get_value(Fields.title)
            
            if not record.has_field(Fields.cycle_id):

                schedule = size_dependent_schedule_cycles(opt, info_dic)
                info_str = str_schedule(schedule, 0, flask_title)
                opt["Logger"].info(info_str)

                new_dic = {
                    "time": schedule[0][0],
                    "CubeTitle": "cycle_0",
                    "suffix": "cycle_0",
                    "trajectory_interval": opt["trajectory_interval"],
                    "reporter_interval": opt["reporter_interval"],
                    "hmr": info_dic["hmr"],
                    "md_engine": info_dic["md_engine"],
                    "instance_type": info_dic["instance_type"],
                    "save_md_stage": True
                }

                record.set_value(Fields.cycle_id, 0)
                record.set_value(Fields.schedule, schedule)
                record.set_value(Fields.cube_parameters_update, new_dic)

                Terminate=False
                time_remain= opt["time"] - new_dic["time"]

            else: 

                current_cycle_id = record.get_value(Fields.cycle_id)
                schedule = record.get_value(Fields.schedule)
                next_cycle_id = current_cycle_id + 1

                # Cycle Termination
                if str(next_cycle_id) not in schedule and not record.has_field(Fields.end_cycle):
                    opt["Logger"].info("Cycle complete....{}".format(flask_title))
                    record.set_value(Fields.cycle_id, next_cycle_id)
                    record.set_value(Fields.end_cycle, True)

                # Cycle Restarting
                elif str(next_cycle_id) not in schedule and record.has_field(Fields.end_cycle):
                    opt["Logger"].info("Restarting Detected....")                

                    schedule = size_dependent_schedule_cycles(opt, info_dic, start_iter_index=current_cycle_id)
                    info_str = str_schedule(schedule, current_cycle_id, flask_title)
                    opt["Logger"].info(info_str)

                    new_dic = {
                        "time": schedule[current_cycle_id][0],
                        "CubeTitle": "cycle_" + str(current_cycle_id),
                        "suffix": "cycle_" + str(current_cycle_id),
                        "trajectory_interval": opt["trajectory_interval"],
                        "reporter_interval": opt["reporter_interval"],
                        "hmr": info_dic["hmr"],
                        "md_engine": info_dic["md_engine"],
                        "instance_type": info_dic["instance_type"],
                        "save_md_stage": True}

                    if record.has_field(Fields.end_cycle):
                        record.delete_field(Fields.end_cycle)

                    record.set_value(Fields.schedule, schedule)
                    record.set_value(Fields.cycle_id, current_cycle_id)
                    record.set_value(Fields.cube_parameters_update, new_dic)

                    Terminate=False
                    time_remain= opt["time"] - new_dic["time"]

                # Cycle forwarding
                else:

                    schedule = size_dependent_schedule_cycles(opt, info_dic, start_iter_index=next_cycle_id)
                    info_str = str_schedule(schedule, next_cycle_id, flask_title)
                    opt["Logger"].info(info_str)

                    new_dic = {
                        "time": schedule[next_cycle_id][0],
                        "CubeTitle": "cycle_" + str(next_cycle_id),
                        "suffix": "cycle_" + str(next_cycle_id),
                        "trajectory_interval": opt["trajectory_interval"],
                        "reporter_interval": opt["reporter_interval"],
                        "hmr": info_dic["hmr"],
                        "md_engine": info_dic["md_engine"],
                        "instance_type": info_dic["instance_type"],
                        "save_md_stage": True
                    }

                    record.set_value(Fields.schedule, schedule)
                    record.set_value(Fields.cycle_id, next_cycle_id)
                    record.set_value(Fields.cube_parameters_update, new_dic)

                    Terminate=False
                    time_remain= opt["time"] - new_dic["time"]                    

            ##############################
            # 3. Run sincle cycle of NPT #
            ##############################

            if Terminate is False: 

                # Update Cube Parameters
                update_cube_parameters_in_place(record, opt)

                # Run the MD simulation
                mdstate = mdrecord.get_stage_state()
                parmed_structure = mdrecord.get_parmed(sync_stage_name="last")

                new_mdstate = md_simulation(mdstate, parmed_structure, opt)

                # Update the system coordinates
                flask.SetCoords(new_mdstate.get_oe_positions())
                mdrecord.set_flask(flask)

                # Trajectory
                if opt["trajectory_interval"] or opt["trajectory_frames"]:
                    trajectory_fn = opt["trj_fn"]
                    if opt["md_engine"] == MDEngines.OpenMM:
                        trajectory_engine = MDEngines.OpenMM
                        traj_ext = '.h5'
                    else:
                        trajectory_engine = MDEngines.Gromacs
                        traj_ext = '.gromacs'
                    # HJ: find extension of the trj file
                    trajectory_orion_ui = opt["system_title"] + "_" + str(opt["system_id"]) + "-" + opt["suffix"]+'_traj' + traj_ext + ".tar.gz"

                else:  # Empty Trajectory
                    trajectory_fn = None
                    trajectory_engine = None
                    trajectory_orion_ui = opt["system_title"]+ "_"+ str(opt["system_id"])+ "-"+ opt["suffix"]+ ".tar.gz" # HJ

                data_fn = opt["out_fn"] + ".tar.gz"

                # Save the cube parameters that are serializable
                info_dic = dict()
                for k, v in dict(vars(self.args)).items():
                    if (
                        isinstance(v, str)
                        or isinstance(v, int)
                        or isinstance(v, float)
                        or isinstance(v, bool)
                    ):
                        info_dic[k] = v

                info_dic["speed_ns_per_day"] = opt["speed_ns_per_day"]
                info_dic["CubeTitle"] = opt["CubeTitle"]

                # Update common values
                for k in info_dic:
                    if k in opt:
                        info_dic[k] = opt[k]

                max_attempt = 5

                for try_count in range(0, max_attempt):
                
                    print("....Attempt {} to add new stage: {}".format(try_count + 1, opt["CubeTitle"]))
                    success = False

                    try: 
                        mdrecord.add_new_stage(
                            opt["CubeTitle"],
                            MDStageTypes.NPT,
                            flask,
                            new_mdstate,
                            data_fn,
                            append=opt["save_md_stage"],
                            log=opt["str_logger"],
                            info=info_dic,
                            trajectory_fn=trajectory_fn,
                            trajectory_engine=trajectory_engine,
                            trajectory_orion_ui=trajectory_orion_ui,
                        )
                        # reset cube_parameters_update
                        new_dic = {
                            "time": time_remain
                        }
                        record.set_value(Fields.cube_parameters_update, new_dic)
                        success=True

                    except Exception as e:
                        if mdrecord.has_stage_name(stg_name=opt["CubeTitle"]):
                            mdrecord.delete_stage_by_name(stg_name=opt["CubeTitle"])

                    if success:
                        print("....Done adding new stage {}".format(opt["CubeTitle"]))
                        break

                if try_count == max_attempt - 1:
                    raise ValueError("Failed to add new stage: {}".format(opt["CubeTitle"]))

                self.cycle_out_port.emit(mdrecord.get_record)

                del mdrecord
            
            else: 
                # Cycle Termination 
                self.success.emit(record)

        except Exception as e:
            msg = "{} Cube exception: {}".format(self.title, str(e))
            self.opt["Logger"].info(msg)
            self.log.error(traceback.format_exc())
            # Write field for Exception Handler
            record.set_value(self.args.log_field, str(msg))
            # Return failed mol
            self.failure.emit(record)

        return
    
  
class ParallelIterMDNptCube(ParallelMixin, IterMDNptCube):
    title = "Parallel " + IterMDNptCube.title
    description = "(Parallel) " + IterMDNptCube.description
    uuid = "4a784594-e887-4748-9362-2f792f29e6b4"

    # Override defaults for some parameters
    parameter_overrides = {
        "gpu_count": {"default": 1},
        "instance_type": {"default": GPU_INSTANCE_TYPES_SUPPORT_FOR_MD},  # Gpu Family selection
        # "instance_type": {"default": ""},  # Gpu Family selection
        "memory_mb": {"default": 14000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1},  # 1 molecule at a time
    }

class IterChecker(RecordPortsMixin, ComputeCube):
    title = "Iteration Checker Cube"
    classification = [["MD Simulations"]]
    tags = ["Gromacs", "OpenMM"]
    description="""
    This Cube checks if the MD cycle is complete. If the MD cycle is complete, 
    it emits the record to the "true" output port and if not, it emits the record 
    to the "false" output port. 
    """
    related = ["BooleanFieldSwitch"]
    uuid = "049765be-de4f-4c58-81a5-5bf2b5b00ff0"

    # for Exception Handler
    log_field = LogFieldParam()

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 32000},
        "disk_space": {"default": 100000},
        "spot_policy": {"default": "Prohibited"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }
    
    cycle_in_port = RecordInputPort("cycle_in_port")
    intake = RecordInputPort("intake")

    true = RecordOutputPort()    
    false = RecordOutputPort()

    recover_out_port = RecordOutputPort("recover_out_port")

    def begin(self):
        self.opt = vars(self.args)
        self.opt["Logger"] = self.log

        return
    
    def process(self, record, port):
        
        try:
            opt = dict(self.opt)
            opt["CubeTitle"] = self.title 
            
            if port == "cycle_in_port":
                self.true.emit(record)
                self.recover_out_port.emit(record)
                
            elif port == "intake":                
                mdrecord = MDDataRecord(record)
                
                # 1. Get topology, mdstate, info_dic and log from the last state
                last_flask = mdrecord.get_stage_topology(stg_name="last")
                last_mdstate = mdrecord.get_stage_state(stg_name="last")
                last_info_dic = mdrecord.get_stage_info(stg_name="last")
                last_stage_log = mdrecord.get_stage_logs(stg_name="last")

                system_title = mdrecord.get_title
                system_id = mdrecord.get_flask_id
                suffix = "comb"

                # 2. set filenames
                cwd = mdrecord.cwd
                out_directory = cwd
                out_fn = check_filename(
                        out_directory
                        + "_"
                        + system_title + "_"
                        + str(system_id)
                        + "-"
                        + suffix
                )
                data_fn = out_fn + ".tar.gz"
                tar_fn = out_fn + "_" + "traj.tar.gz"

                # 3. get trajectory_engine and trajectory_orion_ui
                if last_info_dic["md_engine"] == MDEngines.OpenMM:
                    trajectory_engine = MDEngines.OpenMM
                    traj_ext = '.h5'
                else:
                    trajectory_engine = MDEngines.Gromacs
                    traj_ext = '.gromacs'
                if last_info_dic["trajectory_interval"] or last_info_dic["trajectory_frames"]:
                    trajectory_orion_ui = system_title + "_" + str(system_id) + "-" + suffix +'_traj' + traj_ext + ".tar.gz"
                else:  # Empty Trajectory
                    tar_fn = None
                    trajectory_engine = None
                    trajectory_orion_ui = system_title + "_" + str(system_id) + "-" + suffix + ".tar.gz"

                stage_names = [stg_name for stg_name in mdrecord.get_stages_names if
                                'Production' in stg_name or 'cycle_' in stg_name]
                
                # 4. Combine trajectories into one
                trj_fn_lst = list()
                log_fn_lst = list()
                edr_fn_lst = list()
                tpr_fn_lst = list()
                for idx, stg_name in enumerate(stage_names):
                    trj_fn = mdrecord.get_stage_trajectory(stg_name=stg_name)
                    void, traj_ext = os.path.splitext(trj_fn)
                    traj_dir = os.path.dirname(trj_fn)

                    trj_fn_lst.append(trj_fn)

                    if traj_ext == '.trr':
                        # from the last state
                        gro_fn = glob.glob(os.path.join(traj_dir, '*.gro'))[0]
                        pdb_fn = glob.glob(os.path.join(traj_dir, '*.pdb'))[0]
                        mdp_fn = glob.glob(os.path.join(traj_dir, '*.mdp'))[0] # md param file
                        top_fn = glob.glob(os.path.join(traj_dir, '*.top'))[0] # topology file
                        tpr_fn = glob.glob(os.path.join(traj_dir, '*.tpr'))[0] # input file
                        log_fn = glob.glob(os.path.join(traj_dir, '*.log'))[0] # log file
                        edr_fn = glob.glob(os.path.join(traj_dir, '*.edr'))[0] # energy file

                        log_fn_lst.append(log_fn)
                        edr_fn_lst.append(edr_fn)
                        tpr_fn_lst.append(tpr_fn)


                # 5. Write tar_fn
                if traj_ext == '.h5':
                    omm_trj_fn = os.path.join(out_directory, "trajectory_comb.h5")
                    # combine h5 files into one

                    def combine_h5(h5_fnm_lst, out_h5_fnm):
                        import h5py
                        import numpy as np

                        data_dict = dict()
                        attrs_dict = dict()
                        for trj_fn in h5_fnm_lst:
                            f = h5py.File(trj_fn, 'r')
                            fdata = dict()
                            for key in f:
                                fdata[key] = f[key][...]
                                if key not in attrs_dict:
                                    attrs_dict[key] = dict(f[key].attrs)
                            f.close()
                            data_dict[trj_fn] = fdata

                        data = None
                        for f in data_dict:
                            if not data:
                                data = data_dict[f]
                            else:
                                for key in data_dict[f]:
                                    data[key] = np.append(data[key], data_dict[f][key], axis=0)

                        f = h5py.File(out_h5_fnm, 'w')
                        for key in data:
                            f.create_dataset(key, data[key].shape, dtype=data[key].dtype,
                                            compression='gzip')[...] = data[key]
                            for attr in attrs_dict[key].keys():
                                print(attr, attrs_dict[key][attr])
                                f[key].attrs[attr] = attrs_dict[key][attr]
                        f.close()

                        return
                    
                    combine_h5(trj_fn_lst, omm_trj_fn)

                    with tarfile.open(tar_fn, mode="w") as archive:
                            archive.add(
                                omm_trj_fn,
                                arcname=os.path.basename(omm_trj_fn),
                            )

                    # Check trajectory file
                    with tempfile.TemporaryDirectory() as tmp_dir:
                        try:
                            with tarfile.open(tar_fn) as tar:
                                tar.extractall(path=tmp_dir)
                            tmp_trj_fn = os.path.join(tmp_dir, os.path.basename(omm_trj_fn))
                            mdtraj.load(tmp_trj_fn)
                        except:
                            raise ValueError("Trajectory Check Error")
                        
                elif traj_ext == '.trr':

                    # Gromacs file names
                    outfname = check_filename(system_title + "_" + str(system_id) + "-" + suffix)

                    grm_top_fn = os.path.join(out_directory, outfname + ".top")
                    grm_gro_fn = os.path.join(out_directory, outfname + ".gro")
                    grm_pdb_fn = os.path.join(out_directory, outfname + ".pdb")
                    grm_mdp_fn = os.path.join(out_directory, outfname + ".mdp")
                    grm_trj_fn = os.path.join(out_directory, outfname + ".trr")

                    # combine trr files into one
                    trrs_to_comb = " ".join(str(trr_fn) for trr_fn in trj_fn_lst)
                    os.system('gmx trjcat -f %s -o %s -cat yes' % (trrs_to_comb, grm_trj_fn))

                    with tarfile.open(tar_fn, mode="w") as archive:

                        archive.add(
                            gro_fn,
                            arcname=os.path.basename(grm_gro_fn),
                        )
                        archive.add(
                            pdb_fn,
                            arcname=os.path.basename(grm_pdb_fn),
                        )
                        archive.add(
                            top_fn,
                            arcname=os.path.basename(grm_top_fn),
                        )
                        archive.add(
                            grm_trj_fn,
                            arcname=os.path.basename(grm_trj_fn),
                        )
                        for edr_fn in edr_fn_lst:
                            archive.add(
                                edr_fn,
                                arcname=os.path.basename(edr_fn),
                            )
                        for log_fn in log_fn_lst:
                            archive.add(
                                log_fn,
                                arcname=os.path.basename(log_fn),
                            )
                        for tpr_fn in tpr_fn_lst:
                            archive.add(
                                tpr_fn,
                                arcname=os.path.basename(tpr_fn),
                            )
                        archive.add(
                            mdp_fn, arcname=os.path.basename(grm_mdp_fn)
                        )

                    # Check trajectory file
                    with tempfile.TemporaryDirectory() as tmp_dir:
                        try:
                            with tarfile.open(tar_fn) as tar:
                                tar.extractall(path=tmp_dir)

                            tmp_trj_fn = os.path.join(tmp_dir, os.path.basename(grm_trj_fn))
                            tmp_top_fn = os.path.join(tmp_dir, os.path.basename(grm_pdb_fn))
                            mdtraj.load(tmp_trj_fn, top=tmp_top_fn)
                        except:
                            raise ValueError("Trajectory Check Error")
                            
                # delete the stages 
                for stg_name in stage_names:
                    opt["Logger"].info("[{}] Deleting stage: {}".format(opt["CubeTitle"], stg_name))
                    mdrecord.delete_stage_by_name(stg_name=stg_name)

                # add one new stage
                stage_name="Production"
                opt["Logger"].info("[{}] Adding a combined stage: {}".format(opt["CubeTitle"], stage_name))

                max_attempt = 5

                for try_count in range(0, max_attempt):
                
                    print("....Attempt {} to add a combined stage: {}".format(try_count + 1, opt["CubeTitle"]))
                    
                    try:                     
                        mdrecord.add_new_stage(
                            stage_name=stage_name,
                            stage_type=MDStageTypes.NPT,
                            topology=last_flask,
                            mdstate=last_mdstate,
                            data_fn=data_fn, 
                            append=last_info_dic["save_md_stage"],
                            log=last_stage_log,
                            info=last_info_dic,
                            trajectory_fn=tar_fn,
                            trajectory_engine=trajectory_engine,
                            trajectory_orion_ui=trajectory_orion_ui,
                        )
                        success=True

                    except Exception as e:
                        if mdrecord.has_stage_name(stg_name=stage_name):
                            mdrecord.delete_stage_by_name(stg_name=stage_name)

                    if success:
                        print("....Done adding a combined stage {}".format(stage_name))
                        break

                if try_count == max_attempt - 1:
                    raise ValueError("Failed to add a combined stage: {}".format(stage_name))

                self.false.emit(record)

        except Exception as e:
            msg = "{} Cube exception: {}".format(self.title, str(e))
            self.opt["Logger"].info(msg)
            self.log.error(traceback.format_exc())
            # Write field for Exception Handler
            record.set_value(self.args.log_field, str(msg))
            # Return failed mol
            self.failure.emit(record)

        return
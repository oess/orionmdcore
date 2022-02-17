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
    from openeye import oechem

    from orionplatform.mixins import RecordPortsMixin

    from floe.api import ParallelMixin, parameters, ComputeCube

    from orionplatform.ports import RecordInputPort, RecordOutputPort
except ImportError:
    from orionmdcore import __installation__error__
    raise ImportError(__installation__error__)

import traceback

from orionmdcore.standards.standards import (
    MDStageTypes,
    MDEngines,
    MDStageNames,
    Fields,
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
# TODO REMOVE THIS
from orionmdcore.cubes.md.utils import str_schedule, schedule_cycles
###


class MDMinimizeCube(RecordPortsMixin, ComputeCube):
    title = "Minimization Cube"

    classification = [["MD Simulations"]]
    tags = ["OpenMM", "Gromacs", "Minimization"]

    description = """
    This Cube performs energy minimization on the provided system. The system 
    must have been parametrized by the Force Field cube and the system Parmed 
    object must be present on the input record. In addition, a system identification 
    number must be present on the input record as well. This can be accomplished 
    by using the “ID Setting Cube”. The system minimization is performed by 
    the selected MD engine, currently OpenMM and Gromacs only. Restraints 
    and constraints can be used as well. Currently implicit solvent models can 
    be used in OpenMM only. The cube requires a record as input and produces 
    a new record with the minimized system.
    """

    uuid = "bdfeaabe-f93b-4a14-9754-d6ca0c18a009"

    # Override defaults for some parameters
    parameter_overrides = {
        "gpu_count": {"default": 1},
        "cpu_count": {"default": 16},
        # "instance_type": {"default": "g3.4xlarge"},  # Gpu Family selection
        "instance_type": {"default": ""},  # Gpu Family selection
        "memory_mb": {"default": 14000},
        "spot_policy": {"default": "Prohibited"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1},  # 1 molecule at a time
    }

    steps = parameters.IntegerParameter(
        "steps",
        default=2000,
        help_text="""Number of minimization steps.
                  If 0 the minimization will continue
                  until convergence""",
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
        default=5.0,
        help_text="Restraint weight for xyz atom restraints in kcal/(mol A^2)",
    )

    restraint_to_reference = parameters.BooleanParameter(
        "restraint_to_reference",
        default=True,
        help_text="If True the starting reference system coordinates will be used "
        "to restraint the system",
    )

    freeze = parameters.StringParameter(
        "freeze",
        default="",
        help_text="""Mask selection to freeze atoms along the MD
        simulation. Possible keywords are: ligand, protein, water,
        ions, ca_protein, cofactors. The selection can be refined by
        using logical tokens: not, noh, and, or, diff, around. NOTE:
        Not currently implemented in Gromacs""",
    )

    temperature = parameters.DecimalParameter(
        "temperature", default=300, help_text="Temperature (Kelvin)"
    )

    nonbondedCutoff = parameters.DecimalParameter(
        "nonbondedCutoff",
        default=10,
        help_text="""The non-bonded cutoff in angstroms.
        This is ignored if the non-bonded method is NoCutoff""",
    )

    constraints = parameters.StringParameter(
        "constraints",
        default="Bonds2H",
        choices=["None", "Bonds2H", "Angles2H", "All-Bonds"],
        help_text="""None, Bonds2H, Angles2H, or All-Bonds
        Which type of constraints to add to the system.
        None means no bonds are constrained.
        Bonds2H means bonds with hydrogen are constrained, etc.""",
    )

    implicit_solvent = parameters.StringParameter(
        "implicit_solvent",
        default="None",
        choices=["None", "HCT", "OBC1", "OBC2", "GBn", "GBn2"],
        help_text="Implicit Solvent Model. NOTE:"
        "Not currently implemented in Gromacs",
    )

    center = parameters.BooleanParameter(
        "center",
        default=True,
        description="Center the system to the OpenMM and Gromacs unit cell",
    )

    verbose = parameters.BooleanParameter(
        "verbose", default=True, description="Increase log file verbosity"
    )

    suffix = parameters.StringParameter(
        "suffix", default="min", help_text="Filename suffix for output simulation files"
    )

    hmr = parameters.BooleanParameter(
        "hmr",
        default=False,
        description="On enables Hydrogen Mass Repartitioning. NOTE:"
        "Not currently implemented in Gromacs",
    )

    save_md_stage = parameters.BooleanParameter(
        "save_md_stage",
        default=True,
        help_text="""Save the MD simulation stage. If True the MD,
           simulation data will be appended to the md simulation stages 
           otherwise the last MD stage will be overwritten""",
    )

    md_engine = parameters.StringParameter(
        "md_engine",
        default="OpenMM",
        choices=["OpenMM", "Gromacs"],
        help_text="Select the MD available engine",
    )

    def begin(self):
        self.opt = vars(self.args)
        self.opt["Logger"] = self.log
        self.opt["SimType"] = "min"
        return

    def process(self, record, port):
        try:
            # The copy of the dictionary option as local variable
            # is necessary to avoid filename collisions due to
            # the parallel cube processes
            opt = dict(self.opt)
            opt["CubeTitle"] = self.title

            # Update Cube Parameters
            update_cube_parameters_in_place(record, opt)

            # Logger string
            str_logger = "-" * 32 + " MIN CUBE PARAMETERS " + "-" * 32
            str_logger += "\n{:<25} = {:<10}".format("Cube Title", opt["CubeTitle"])

            for k, v in sorted(self.parameters().items()):
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
                "[{}] MINIMIZING Flask: {}".format(opt["CubeTitle"], system_title)
            )

            # Extract the Parmed structure and synchronize it with the last MD stage state
            parmed_structure = mdrecord.get_parmed(sync_stage_name="last")

            # Run the MD simulation
            new_mdstate = md_simulation(mdstate, parmed_structure, opt)

            # Update the flask coordinates
            flask.SetCoords(new_mdstate.get_oe_positions())
            mdrecord.set_flask(flask)

            data_fn = check_filename(
                os.path.basename(mdrecord.cwd)
                + "_"
                + opt["system_title"]
                + "_"
                + str(opt["system_id"])
                + "-"
                + opt["suffix"]
                + ".tar.gz"
            )

            # Save the cube parameters tha are serializable
            info_dic = dict()
            for k, v in dict(vars(self.args)).items():
                if (
                    isinstance(v, str)
                    or isinstance(v, int)
                    or isinstance(v, float)
                    or isinstance(v, bool)
                ):
                    info_dic[k] = v

            info_dic["CubeTitle"] = opt["CubeTitle"]

            # Update common values that could have been updated
            for k in info_dic:
                if k in opt:
                    info_dic[k] = opt[k]

            if not mdrecord.add_new_stage(
                opt["CubeTitle"],
                MDStageTypes.MINIMIZATION,
                flask,
                new_mdstate,
                data_fn,
                append=opt["save_md_stage"],
                log=opt["str_logger"],
                info=info_dic,
            ):

                raise ValueError("Problems adding the new Minimization Stage")

            self.success.emit(mdrecord.get_record)

            del mdrecord

        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt["Logger"].info("Exception {} {}".format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)

        return


class MDNvtCube(RecordPortsMixin, ComputeCube):
    title = "NVT Cube"

    classification = [["MD Simulations"]]
    tags = ["Gromacs", "OpenMM", "NVT"]

    description = """
    This Cube performs MD simulation in the NVT ensemble on the provided system. 
    The system must have been parametrized by the Force Field cube and the system Parmed 
    object must be present on the input record. In addition, a system identification 
    number must be present on the input record as well. This can be accomplished 
    by using the “ID Setting Cube”. The NVT MD simulation is performed by the selected 
    MD engine, currently OpenMM and Gromacs only. Restraints and constraints can be 
    used as well. Currently implicit solvent models can be used in OpenMM only. 
    The cube requires a record as input and produces a new record with the time evolved 
    system. The total sampling time can be set by using the “time” cube parameter 
    while the trajectory snapshots can be set by using the “trajectory_interval” cube 
    parameter.
    """

    uuid = "94962b93-bf32-4b5e-b324-dbe8b9350266"

    # Override defaults for some parameters
    parameter_overrides = {
        "gpu_count": {"default": 1},
        "cpu_count": {"default": 16},
        # "instance_type": {"default": "g3.4xlarge"},  # Gpu Family selection
        "instance_type": {"default": ""},  # Gpu Family selection
        "memory_mb": {"default": 14000},
        "spot_policy": {"default": "Prohibited"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1},  # 1 molecule at a time
    }

    temperature = parameters.DecimalParameter(
        "temperature", default=300.0, help_text="Temperature (Kelvin)"
    )

    time = parameters.DecimalParameter(
        "time", default=0.01, help_text="NVT simulation time in nanoseconds"
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
        default=True,
        help_text="If True the starting reference system coordinates will be used "
        "to restraint the system",
    )

    freeze = parameters.StringParameter(
        "freeze",
        default="",
        help_text="""Mask selection to freeze atoms along the MD
        simulation. Possible keywords are: ligand, protein, water,
        ions, ca_protein, cofactors. The selection can be refined by
        using logical tokens: not, noh, and, or, diff, around. NOTE:
        Not currently implemented in Gromacs""",
    )

    nonbondedCutoff = parameters.DecimalParameter(
        "nonbondedCutoff",
        default=10,
        help_text="""The non-bonded cutoff in angstroms.
        This is ignored if non-bonded method is NoCutoff.""",
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
        default="None",
        choices=["None", "HCT", "OBC1", "OBC2", "GBn", "GBn2"],
        help_text="Implicit Solvent Model. NOTE:"
        "Not currently implemented in Gromacs",
    )

    trajectory_interval = parameters.DecimalParameter(
        "trajectory_interval",
        default=0.0,
        help_text="""Time interval for trajectory snapshots in ns. 
        If 0 the trajectory file will not be generated""",
    )

    reporter_interval = parameters.DecimalParameter(
        "reporter_interval",
        default=0.0,
        help_text="""Time interval for reporting data in ns. 
        If 0 the reporter file will not be generated""",
    )

    trajectory_frames = parameters.IntegerParameter(
        "trajectory_frames",
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
        "suffix", default="nvt", help_text="Filename suffix for output simulation files"
    )

    center = parameters.BooleanParameter(
        "center", default=False, help_text="Center the system to the OpenMM unit cell"
    )

    verbose = parameters.BooleanParameter(
        "verbose", default=True, help_text="Increase log file verbosity"
    )

    hmr = parameters.BooleanParameter(
        "hmr",
        default=False,
        help_text="On enables Hydrogen Mass Repartitioning. NOTE: "
        "Not currently implemented in Gromacs",
    )

    save_md_stage = parameters.BooleanParameter(
        "save_md_stage",
        default=False,
        help_text="""Save the MD simulation stage. If True the MD,
           simulation data will be appended to the md simulation stages 
           otherwise the last MD stage will be overwritten""",
    )

    md_engine = parameters.StringParameter(
        "md_engine",
        default="OpenMM",
        choices=["OpenMM", "Gromacs"],
        help_text="Select the MD available engine",
    )

    def begin(self):
        self.opt = vars(self.args)
        self.opt["Logger"] = self.log
        self.opt["SimType"] = "nvt"

        return

    def process(self, record, port):
        try:
            # The copy of the dictionary option as local variable
            # is necessary to avoid filename collisions due to
            # the parallel cube processes
            opt = dict(self.opt)
            opt["CubeTitle"] = self.title

            # Update Cube Parameters
            update_cube_parameters_in_place(record, opt)

            # Logger string
            str_logger = "-" * 32 + " NVT CUBE PARAMETERS " + "-" * 32
            str_logger += "\n{:<25} = {:<10}".format("Cube Title", opt["CubeTitle"])

            for k, v in sorted(self.parameters().items()):
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
                "[{}] START NVT SIMULATION: {}".format(opt["CubeTitle"], system_title)
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

            # Run the MD simulation
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

            # Save the cube parameters tha are serializable
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

            # Update common values that could have been updated
            for k in info_dic:
                if k in opt:
                    info_dic[k] = opt[k]

            if not mdrecord.add_new_stage(
                opt["CubeTitle"],
                MDStageTypes.NVT,
                flask,
                new_mdstate,
                data_fn,
                append=opt["save_md_stage"],
                log=opt["str_logger"],
                info=info_dic,
                trajectory_fn=trajectory_fn,
                trajectory_engine=trajectory_engine,
                trajectory_orion_ui=trajectory_orion_ui,
            ):

                raise ValueError("Problems adding in the new NVT Stage")

            self.success.emit(mdrecord.get_record)

            del mdrecord

        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt["Logger"].info("Exception {} {}".format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)

        return


class MDNptCube(RecordPortsMixin, ComputeCube):
    title = "NPT Cube"

    classification = [["MD Simulations"]]
    tags = ["Gromacs", "OpenMM", "NPT"]

    description = """
    This Cube performs MD simulation in the NPT ensemble on the provided system. 
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
    """

    uuid = "602d397b-d8a5-4388-a94a-ac3a54ff3bad"

    # Override defaults for some parameters
    parameter_overrides = {
        "gpu_count": {"default": 1},
        "cpu_count": {"default": 16},
        # "instance_type": {"default": "g3.4xlarge"},  # Gpu Family selection
        "instance_type": {"default": ""},  # Gpu Family selection
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
        "time", default=0.01, help_text="NPT simulation time in nanoseconds"
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
        default="None",
        choices=["None", "HCT", "OBC1", "OBC2", "GBn", "GBn2"],
        help_text="Implicit Solvent Model. Not Currently implemented in Gromacs",
    )

    trajectory_interval = parameters.DecimalParameter(
        "trajectory_interval",
        default=0.0,
        help_text="""Time interval for trajectory snapshots in ns. 
        If 0 the trajectory file will not be generated""",
    )

    reporter_interval = parameters.DecimalParameter(
        "reporter_interval",
        default=0.0,
        help_text="""Time interval for reporting data in ns. 
        If 0 the reporter file will not be generated""",
    )

    trajectory_frames = parameters.IntegerParameter(
        "trajectory_frames",
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
        default=True,
        help_text="On enables Hydrogen Mass Repartitioning. Not currently implemented in Gromacs",
    )

    save_md_stage = parameters.BooleanParameter(
        "save_md_stage",
        default=False,
        help_text="""Save the MD simulation stage. If True the MD,
           simulation data will be appended to the md simulation stages 
           otherwise the last MD stage will be overwritten""",
    )

    md_engine = parameters.StringParameter(
        "md_engine",
        default="OpenMM",
        choices=["OpenMM", "Gromacs"],
        help_text="Select the MD available engine",
    )

    def begin(self):
        self.opt = vars(self.args)
        self.opt["Logger"] = self.log
        self.opt["SimType"] = "npt"

        return

    def process(self, record, port):
        try:
            # The copy of the dictionary option as local variable
            # is necessary to avoid filename collisions due to
            # the parallel cube processes
            opt = dict(self.opt)
            opt["CubeTitle"] = self.title

            # Update Cube Parameters
            update_cube_parameters_in_place(record, opt)

            # Logger string
            str_logger = "-" * 32 + " NPT CUBE PARAMETERS " + "-" * 32
            str_logger += "\n{:<25} = {:<10}".format("Cube Title", opt["CubeTitle"])

            for k, v in sorted(self.parameters().items()):
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
                "[{}] START NPT SIMULATION: {}".format(opt["CubeTitle"], system_title)
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

            # Run the MD simulation
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

            if not mdrecord.add_new_stage(
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
            ):

                raise ValueError("Problems adding in the new NPT Stage")

            self.success.emit(mdrecord.get_record)

            del mdrecord

        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt["Logger"].info("Exception {} {}".format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)

        return


class MDProxyCube(RecordPortsMixin, ComputeCube):
    title = "MD ProxyCube"

    classification = [["MD Simulations"]]
    tags = ["Gromacs", "OpenMM"]

    description = """
    The MD Proxy cube runs chucks of MD simulations in cycles till the desired
    total md running time has been reached. 

    Input:
    -------
    Data Record Streaming - Streams of MD ready records 

    Output:
    -------
    Data Record Streaming - Streams of MD ready records with the long trajectory md
    stages in it

    Data Record Streaming - Streams of Recovery records that can be used in the event 
    of Orion issues
    """

    uuid = "1c0882c5-be6c-4405-8fba-ef78cdbc2429"

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 14000},
        "spot_policy": {"default": "Prohibited"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1},  # 1 molecule at a time
    }

    time = parameters.DecimalParameter(
        "time", default=5.0, help_text="Total simulation time in nanoseconds"
    )

    trajectory_interval = parameters.DecimalParameter(
        "trajectory_interval",
        default=0.10,
        help_text="""Time interval per cycle for trajectory snapshots in ns. 
        If 0 the trajectory file will not be generated""",
    )

    reporter_interval = parameters.DecimalParameter(
        "reporter_interval",
        default=0.01,
        help_text="""Time interval for reporting data in ns. 
        If 0 the reporter file will not be generated""",
    )

    cube_max_run_time = parameters.DecimalParameter(
        "cube_max_run_time", default=10, help_text="Max Cube Running Time in hrs. Max 11hrs"
    )

    cycle_out_port = RecordOutputPort(
        "cycle_out_port",
    )

    recover_out_port = RecordOutputPort(
        "recover_out_port",
    )

    def begin(self):
        self.opt = vars(self.args)
        self.opt["Logger"] = self.log

        return

    def process(self, record, port):
        try:
            opt = dict(self.opt)
            opt["CubeTitle"] = self.title

            if not 0 < opt['cube_max_run_time'] < 11:
                raise ValueError("The Max allowed MD running time per iteration is limited to 11 hrs. "
                                 "Selected time: {} hrs".format(opt['cube_max_run_time']))

            md_record = MDDataRecord(record)

            info_dic = md_record.get_stage_info(stg_name="last")

            if not record.has_field(Fields.title):
                raise ValueError("Field title has not been found on the record")

            flask_title = record.get_value(Fields.title)

            if not record.has_field(Fields.cycle_id):

                if "speed_ns_per_day" not in info_dic:
                    raise ValueError("Last MD stage does not have speed information")

                # Remove the init stage
                md_record.delete_stage_by_name(stg_name="last")

                schedule = schedule_cycles(opt, info_dic)
                info_str = str_schedule(schedule, 0, flask_title)
                opt["Logger"].info(info_str)

                record.set_value(Fields.cycle_id, 0)
                record.set_value(Fields.schedule, schedule)

                # TODO UPDATE ALL THE CUBE PARAMS?
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

                record.set_value(Fields.cube_parameters_update, new_dic)
                opt["Logger"].info("Forwarding to next cycle....{}".format(flask_title))
                self.cycle_out_port.emit(record)

            else:

                current_cycle_id = record.get_value(Fields.cycle_id)
                schedule = record.get_value(Fields.schedule)
                next_cycle_id = current_cycle_id + 1

                # Cycle Termination
                if str(next_cycle_id) not in schedule and not record.has_field(Fields.end_cycle):
                    opt["Logger"].info("Cycle complete....{}".format(flask_title))
                    record.set_value(Fields.cycle_id, next_cycle_id)
                    self.recover_out_port.emit(record)
                    record.set_value(Fields.end_cycle, True)
                    self.success.emit(record)

                # Cycle Restarting
                elif str(next_cycle_id) not in schedule and record.has_field(Fields.end_cycle):
                    opt["Logger"].info("Restarting Detected....")

                    # Remove the init stage
                    md_record.delete_stage_by_name(stg_name="last")

                    if "speed_ns_per_day" not in info_dic:
                        raise ValueError("Last MD stage does not have speed information")

                    schedule = schedule_cycles(opt, info_dic, start_iter_index=current_cycle_id)

                    # TODO UPDATE ALL THE CUBE PARAMS?
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

                    record.set_value(Fields.schedule, schedule)
                    record.set_value(Fields.cycle_id, current_cycle_id)
                    record.set_value(Fields.cube_parameters_update, new_dic)

                    info_str = str_schedule(schedule, current_cycle_id, flask_title)
                    opt["Logger"].info(info_str)

                    if record.has_field(Fields.end_cycle):
                        record.delete_field(Fields.end_cycle)

                    self.cycle_out_port.emit(record)

                # Cycle forwarding
                else:

                    # TODO UPDATE ALL THE CUBE PARAMS?
                    new_dic = {
                        "time": schedule[str(next_cycle_id)][0],
                        "CubeTitle": "cycle_" + str(next_cycle_id),
                        "suffix": "cycle_" + str(next_cycle_id),
                        "trajectory_interval": opt["trajectory_interval"],
                        "reporter_interval": opt["reporter_interval"],
                        "hmr": info_dic["hmr"],
                        "md_engine": info_dic["md_engine"],
                        "instance_type": info_dic["instance_type"],
                        "save_md_stage": True
                    }

                    record.set_value(Fields.cycle_id, next_cycle_id)
                    record.set_value(Fields.cube_parameters_update, new_dic)

                    info_str = str_schedule(schedule, next_cycle_id, flask_title)
                    opt["Logger"].info(info_str)
                    opt["Logger"].info("Forwarding to next cycle....{}".format(flask_title))

                    self.cycle_out_port.emit(record)
                    self.recover_out_port.emit(record)

        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt["Logger"].info("Exception {} {}".format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)

        return


class MDRecoveryRestartProxyCube(RecordPortsMixin, ComputeCube):
    title = "MD Recover Restart for ProxyCube"

    classification = [["MD Simulations"]]
    tags = ["Gromacs", "OpenMM"]

    description = """
    This cube is able to restart or recover a previously run long MD simulation.
    This cube is used as support of the MD proxy cube.

    Input:
    -------
    Data Record Streaming - Streams of MD ready records 

    Output:
    -------
    Data Record Streaming - Streams of MD ready records 
    """

    uuid = "5be8a53a-6de0-492c-9564-3957d88084f2"

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 14000},
        "spot_policy": {"default": "Prohibited"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1},  # 1 molecule at a time
    }

    def begin(self):
        self.opt = vars(self.args)
        self.opt["Logger"] = self.log
        self.max_cycles_dic = dict()

        return

    def process(self, record, port):
        try:

            md_record = MDDataRecord(record)

            if not md_record.has_stages:
                if record.has_field(Fields.oetrajconf_rec):
                    oetraj_conf_records = record.get_value(Fields.oetrajconf_rec)

                    last_conf_rec = oetraj_conf_records[-1]
                    md_last_conf_rec = MDDataRecord(last_conf_rec)

                    if md_last_conf_rec.has_stages:
                        record = last_conf_rec
                        md_record = md_last_conf_rec
                        self.opt["Logger"].info("Starting From MD Analysis Record")
                    else:
                        raise ValueError("The record is missing the MD Stages")
                else:
                    raise ValueError("The record is missing the Conf Records")

            # md_components = md_record.get_md_components
            #
            # # if isinstance(md_components, MDComponents):
            # #     raise ValueError("\nDeprecation Error: The MDComponents object seems to be generated "
            # #                      "by using the old MD API. You should use the MD API Converter floe "
            # #                      "to generate datasets compatible with the new MD API\n")

            if record.has_field(Fields.flaskid):
                flask_id = record.get_value(Fields.flaskid)

                # Start
                if not record.has_field(Fields.cycle_id):
                    self.success.emit(record)
                    self.opt["Logger"].info("Starting From MD Record")
                else:
                    cycle_id = record.get_value(Fields.cycle_id)

                    if not record.has_field(Fields.schedule):
                        raise ValueError("The Schedule Field has not been found")

                    schedule = record.get_value(Fields.schedule)

                    # Restarting
                    if str(cycle_id) not in schedule and record.has_field(Fields.end_cycle):
                        # Remove update cube parameters
                        if record.has_field(Fields.cube_parameters_update):
                            record.delete_field(Fields.cube_parameters_update)
                        self.opt["Logger"].info("\nAttempting Restarting....\n")
                        self.success.emit(record)
                        return

                    # Recovering
                    if flask_id not in self.max_cycles_dic:
                        self.max_cycles_dic[flask_id] = record
                    else:

                        max_cycle_id = self.max_cycles_dic[flask_id].get_value(Fields.cycle_id)
                        if cycle_id > max_cycle_id:
                            self.max_cycles_dic[flask_id] = record
            else:
                raise ValueError("Cannot find the flask id field")

        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt["Logger"].info("Exception {} {}".format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)

        return

    def end(self):
        try:
            if self.max_cycles_dic:

                for rec in self.max_cycles_dic.values():

                    if not rec.has_field(Fields.title):
                        raise ValueError("The Title Field has not been found")

                    schedule = rec.get_value(Fields.schedule)
                    title = rec.get_value(Fields.title)
                    next_it = rec.get_value(Fields.cycle_id)

                    # Recovering from a finished run...Do not do anything
                    if str(next_it) not in schedule and not rec.has_field(Fields.end_cycle):
                        self.opt["Logger"].info(
                            "\nAttempting Restarting from completed Recovery record...{}".format(title))
                        return

                    info_str = str_schedule(schedule, next_it, title)
                    self.opt["Logger"].warn("\nAttempting Recovering {}....\n".format(title))
                    self.opt["Logger"].info(info_str)
                    self.success.emit(rec)

        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt["Logger"].info("Exception {}".format(str(e)))
            self.log.error(traceback.format_exc())


class ParallelMDMinimizeCube(ParallelMixin, MDMinimizeCube):
    title = "Parallel " + MDMinimizeCube.title
    description = "(Parallel) " + MDMinimizeCube.description
    uuid = "24ed6b12-a426-4e0d-bbac-99b9e34cca5c"

    parameter_overrides = {
        "gpu_count": {"default": 1},
        # "instance_type": {"default": "g3.4xlarge"},  # Gpu Family selection
        "instance_type": {"default": ""},  # Gpu Family selection
        "memory_mb": {"default": 14000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1},  # 1 molecule at a time
    }


class ParallelMDNvtCube(ParallelMixin, MDNvtCube):
    title = "Parallel " + MDNvtCube.title
    description = "(Parallel) " + MDNvtCube.description
    uuid = "1cff32be-9b10-4070-a4d2-f7370cc8be96"

    # Override defaults for some parameters
    parameter_overrides = {
        "gpu_count": {"default": 1},
        # "instance_type": {"default": "g3.4xlarge"},  # Gpu Family selection
        "instance_type": {"default": ""},  # Gpu Family selection
        "memory_mb": {"default": 14000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1},  # 1 molecule at a time
    }


class ParallelMDNptCube(ParallelMixin, MDNptCube):
    title = "Parallel " + MDNptCube.title
    description = "(Parallel) " + MDNptCube.description
    uuid = "94728422-e840-49ba-9006-f6170dad54ba"

    # Override defaults for some parameters
    parameter_overrides = {
        "gpu_count": {"default": 1},
        # "instance_type": {"default": "g3.4xlarge"},  # Gpu Family selection
        "instance_type": {"default": ""},  # Gpu Family selection
        "memory_mb": {"default": 14000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1},  # 1 molecule at a time
    }

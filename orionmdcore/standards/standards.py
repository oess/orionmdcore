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
    from floe.api.orion import in_orion
except ImportError:
    from orionmdcore import __installation__error__
    raise ImportError(__installation__error__)

from orionmdcore.standards.utils import (
    ParmedData,
    MDStateData,
    DesignUnit,
    MDComponentData,
)

from datarecord import OEPrimaryMolField

from datarecord import Types, Meta, OEFieldMeta, OEField


# ------------ Stage Standard Names ------------- #


class MDStageTypes:
    SETUP = "SETUP"
    MINIMIZATION = "MINIMIZATION"
    NVT = "NVT"
    NPT = "NPT"
    FEC = "FEC"


class MDStageNames:
    ForceField = "Flask Parametrization"
    Minimization = "Flask Minimization"
    WarmUp = "WarmUp"
    EquilibrationI = "EquilibrationI"
    EquilibrationII = "EquilibrationII"
    EquilibrationIII = "EquilibrationIII"
    EquilibrationIV = "EquilibrationIV"
    Production = "Production"


# ------------ MD Engines ------------- #


class MDEngines:
    OpenMM = "OpenMM"
    Gromacs = "Gromacs"
    all = [OpenMM, Gromacs]


# ---------------- File  Name Standards -------------- #


class MDFileNames:
    topology = "topology.oeb"
    state = "state.pickle"
    trajectory = "trajectory.tar.gz"
    trajectory_conformers = "trajectory_confs.oeb"
    mddata = "data.tar.gz"


# ---------------- Collection  Name Standards -------------- #


class CollectionsNames:
    none = ""
    md = "MD_OPLMD"
    nes = "NES_OPLMD"


# Orion Hidden meta data options
_metaHidden = OEFieldMeta(options=[Meta.Display.Hidden])
_metaIDHidden = OEFieldMeta(options=[Meta.Display.Hidden], attributes=[[Meta.Source.ID, "ID"]])
_metaProtHidden = OEFieldMeta(options=[Meta.Hints.Chem.Protein, Meta.Display.Hidden])


# ---------------- Field Standards -------------- #
class Fields:

    # The LigInitialRecord Field is for the initial ligand record read in at the start
    ligInit_rec = OEField("LigInitial", Types.Record, meta=_metaHidden)

    # The Title field is a string name for the flask which used to compose file names
    title = OEField("Title_OPLMD", Types.String, meta=_metaIDHidden)

    # The flaskid field is a unique integer for each flask (final system for simulation)
    flaskid = OEField("FlaskID_OPLMD", Types.Int, meta=_metaIDHidden)

    # The ligid field is a unique integer used to keep track of the ligand input order
    ligid = OEField("LigID_OPLMD", Types.Int, meta=_metaIDHidden)

    # The ConfID field is used to identify a particular conformer
    confid = OEField("ConfID_OPLMD", Types.Int, meta=_metaIDHidden)

    # The Ligand field should be used to save in a record a ligand as an OEMolecule
    ligand = OEField(
        "Ligand_OPLMD",
        Types.Chem.Mol,
        meta=OEFieldMeta(options=[Meta.Hints.Chem.Ligand, Meta.Display.Hidden]),
    )

    # The ligand name
    ligand_name = OEField("Ligand_name_OPLMD", Types.String, meta=_metaHidden)

    # The protein field should be used to save in a record a Protein as an OEMolecule
    protein = OEField("Protein_OPLMD", Types.Chem.Mol, meta=_metaProtHidden)

    # The protein name
    protein_name = OEField("Protein_name_OPLMD", Types.String, meta=_metaHidden)

    # The super-molecule for the entire flask (ie the final system for simulation)
    flask = OEField("Flask_OPLMD", Types.Chem.Mol, meta=_metaHidden)

    # Primary Molecule
    primary_molecule = OEPrimaryMolField()

    # Cube Parameters update record
    cube_parameters_update = OEField(
        "Cube_Parameters_Update_OPLMD", Types.JSONObject, meta=_metaHidden
    )

    # Parmed Structure, Trajectory, MDData and Protein trajectory conformers Fields
    if in_orion():
        pmd_structure = OEField("Structure_Parmed_OPLMD", Types.Int, meta=_metaHidden)
        trajectory = OEField("Trajectory_OPLMD", Types.Int, meta=_metaHidden)
        mddata = OEField("MDData_OPLMD", Types.Int, meta=_metaHidden)

        protein_traj_confs = OEField("ProtTraj_OPLMD", Types.Int, meta=_metaHidden)
        ligand_traj_confs = OEField("LigandTraj_OPLMD", Types.Int, meta=_metaHidden)
        water_traj_confs = OEField("WaterTraj_OPLMD", Types.Int, meta=_metaHidden)

        extra_data_tar = OEField("ExtraData_OPLMD", Types.Int, meta=_metaHidden)
    else:
        pmd_structure = OEField("Structure_Parmed_OPLMD", ParmedData, meta=_metaHidden)
        trajectory = OEField("Trajectory_OPLMD", Types.String, meta=_metaHidden)
        mddata = OEField("MDData_OPLMD", Types.String, meta=_metaHidden)

        protein_traj_confs = OEField("ProtTraj_OPLMD", Types.Chem.Mol, meta=_metaHidden)
        ligand_traj_confs = OEField("LigandTraj_OPLMD", Types.Chem.Mol, meta=_metaHidden)
        water_traj_confs = OEField("WaterTraj_OPLMD", Types.Chem.Mol, meta=_metaHidden)

        extra_data_tar = OEField("ExtraData_OPLMD", Types.String, meta=_metaHidden)

    # The Stage Name
    stage_name = OEField("Stage_name_OPLMD", Types.String)

    # The Stage Type
    stage_type = OEField("Stage_type_OPLMD", Types.String)

    # Topology Field
    topology = OEField(
        "Topology_OPLMD",
        Types.Chem.Mol,
        meta=OEFieldMeta().set_option(Meta.Hints.Chem.PrimaryMol),
    )

    # Log Info
    log_data = OEField("Log_data_OPLMD", Types.String)

    # MD State
    md_state = OEField("MDState_OPLMD", MDStateData)

    # MD Stage Info dictionary
    stage_info = OEField("Stage_info_OPLMD", Types.JSONObject)

    # Design Unit Field
    design_unit = OEField("Design_Unit_OPLMD", DesignUnit)

    # Design Unit Field from Spruce
    # design_unit_from_spruce = OEField('du_single', Types.Blob)
    design_unit_from_spruce = OEField("designunit", Types.Chem.DesignUnit)

    # MD Components
    md_components = OEField("MDComponents_OPLMD", MDComponentData)

    # Collection is used to offload data from the record which must be < 100Mb
    # collection = OEField("Collection_ID_OPLMD", Types.Int, meta=_metaHidden)

    collections = OEField("Collections_ID_OPLMD", Types.JSONObject, meta=_metaHidden)

    cycle_id = OEField("Cycle_ID_OPLMD", Types.Int, meta=_metaHidden)
    end_cycle = OEField("Cycle_END_OPLMD", Types.Bool, meta=_metaHidden)
    schedule = OEField("Schedule_IDS_OPLMD", Types.JSONObject, meta=_metaHidden)

    # Stage list Field
    md_stages = OEField("MDStages_OPLMD", Types.RecordVec, meta=_metaHidden)

    # Bound or Unbound flask type used for Bound or Unbound protein ligand calculations
    bound_unbound_type = OEField("Flask_type", Types.String, meta=_metaHidden)

    floe_report = OEField("Floe_report_OPLMD", Types.String, meta=_metaHidden)

    floe_report_svg_lig_depiction = OEField(
        "Floe_report_lig_svg_OPLMD",
        Types.String,
        meta=OEFieldMeta().set_option(Meta.Hints.Image_SVG),
    )

    floe_report_sort_string = OEField(
        "Floe_report_sort_str_OPLMD", Types.String, meta=_metaHidden
    )

    floe_report_sort_float = OEField(
        "Floe_report_sort_float_OPLMD", Types.Float, meta=_metaHidden
    )

    floe_report_label = OEField(
        "Floe_report_label_OPLMD", Types.String, meta=_metaHidden
    )

    floe_report_URL = OEField(
        "Floe_report_URL_OPLMD",
        Types.String,
        meta=OEFieldMeta(options=[Meta.Hints.URL]),
    )

    floe_report_collection_id = OEField(
        "Floe_report_ID_OPLMD", Types.Int, meta=_metaHidden
    )


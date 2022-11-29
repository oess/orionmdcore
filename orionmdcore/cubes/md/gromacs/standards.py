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
    from datarecord import Types, OEField, OEFieldMeta, Meta
    from orionclient.session import in_orion
except ImportError:
    from orionmdcore import __installation__error__
    raise ImportError(__installation__error__)

# Orion Hidden meta data options
_metaHidden = OEFieldMeta(options=[Meta.Display.Hidden])
_metaIDHidden = OEFieldMeta(options=[Meta.Display.Hidden], attributes=[[Meta.Source.ID, "ID"]])


class Fields:
    # Current number of MD steps
    current_iteration_field = OEField("Current_Iterations_OMD", Types.Int)

    # Total number of MD steps
    md_nsteps_field = OEField("MD_nsteps_OMD", Types.Int)

    # Current number of cycles
    cycle_id = OEField("Cycle_ID_OMD", Types.Int)

    # Tpr binary file
    tpr_field = OEField("TPR_bytes_OMD", Types.Blob, meta=_metaHidden)

    # Prefix name field
    prefix_name_field = OEField("Prefix_OPLMD", Types.String)

    if in_orion():
        trajectory = OEField("GMXTrajectory_OMD", Types.Int, meta=_metaHidden)
        gmx_restart = OEField("GMXRestart_OMD", Types.Int, meta=_metaHidden)
    else:
        trajectory = OEField("GMXTrajectory_OMD", Types.String, meta=_metaHidden)
        gmx_restart = OEField("GMXRestart_OMD", Types.String, meta=_metaHidden)


class Gromacs:
    gmx_tpr_prefix = "gmx_tpr"

    gmx_mdp_prefix = "gmx_mdp"

    gmx_restart_prefix = "gmx_restart"

    gmx_run_prefix = "gmx_run"

    gmx_traj_dir_prefix = "gmx_traj"

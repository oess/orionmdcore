#!/usr/bin/env python

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

from os import path

from floe.api import (WorkFloe)

from orionplatform.cubes import DatasetReaderCube, DatasetWriterCube

from orionmdcore.cubes.flask import (CollectionSetting,
                                     ParallelRecordSizeCheck,
                                     IDSettingCube,
                                     ParallelSolvationCube)

from orionmdcore.cubes.ligprep import (LigandSetting,
                                       ParallelLigandChargeCube)

from snowball import ExceptHandlerCube, DatasetReaderOptCube

from orionmdcore.cubes.flask import MDComponentCube

from orionmdcore.cubes.complexprep import ComplexPrepCube

from orionmdcore.cubes.forcefield import ParallelForceFieldCube

from orionmdcore.cubes.md import (ParallelMDMinimizeCube,
                                  ParallelMDNvtCube,
                                  ParallelMDNptCube)

from floe.api import ParallelCubeGroup


floe_title = 'Bound Protein-Ligand MD'
tags_for_floe = ['MDPrep', 'MDRun']
#
tag_str = ''.join(' [{}]'.format(tag) for tag in tags_for_floe)
job = WorkFloe(floe_title, title=floe_title+tag_str)
job.classification = [tags_for_floe]
job.tags = tags_for_floe

job.description = open(path.join(path.dirname(__file__), 'ProteinLigandMD_desc.rst'), 'r').read()

job.uuid = "ae561d76-a2b6-4d89-b621-b979f1930b40"


# This Cube is necessary for the correct work of collection and shard
coll_open = CollectionSetting("OpenCollection", title="Open Collection")
coll_open.set_parameters(open=True)
coll_open.set_parameters(write_new_collection='MD_OPLMD')

# This Cube is necessary for the correct working of collection and shard
coll_close = CollectionSetting("CloseCollection", title="Close Collection")
coll_close.set_parameters(open=False)

check_rec = ParallelRecordSizeCheck("Record Check Success")

exceptions = ExceptHandlerCube(floe_report_name="Analyze Floe Failure Report")

ofs = DatasetWriterCube('ofs', title='MD Out')
ofs.promote_parameter("data_out", promoted_name="out",
                      title="MD Out", description="MD Dataset out", order=2)

fail = DatasetWriterCube('fail', title='Failures')
fail.promote_parameter("data_out", promoted_name="fail", title="Failures",
                       description="MD Dataset Failures out", order=3)

# Ligand setting
iligs = DatasetReaderCube("LigandReader", title="Ligand Reader")
iligs.promote_parameter("data_in", promoted_name="ligands", title="Ligand Input Dataset",
                        description="Ligand Dataset", order=1)

ligset = LigandSetting("LigandSetting", title="Ligand Setting")
ligset.promote_parameter('max_md_runs', promoted_name='max_md_runs',
                         default=500,
                         description='The maximum allowed number of md runs')
ligset.promote_parameter('n_md_starts', promoted_name='n_md_starts',
                         default=1,
                         description='The number of md starts for each ligand/conformer')
ligset.set_parameters(lig_res_name='LIG')

chargelig = ParallelLigandChargeCube("LigCharge", title="Ligand Charge")
chargelig.promote_parameter('charge_ligands', promoted_name='charge_ligands',
                            description="Charge the ligand or not", default=True)

ligid = IDSettingCube("Ligand Ids")

# Protein Reading cube. The protein prefix parameter is used to select a name for the
# output system files
iprot = DatasetReaderOptCube("ProteinReader", title="Protein Reader")
iprot.promote_parameter("data_in", promoted_name="protein", title='Protein Input Dataset',
                        description="Protein Dataset", order=0)

# Protein Setting
mdcomp = MDComponentCube("MD Components", title="MD Components")
mdcomp.promote_parameter("flask_title", promoted_name="flask_title", default="")

# Complex cube used to assemble the ligands and the solvated protein
complx = ComplexPrepCube("Complex", title="Complex Preparation")

# The solvation cube is used to solvate the system and define the ionic strength of the solution
solvate = ParallelSolvationCube("Solvation", title="Solvation")

# Force Field Application
ff = ParallelForceFieldCube("ForceField", title="Apply Force Field")
ff.promote_parameter('protein_forcefield', promoted_name='protein_ff', default='Amber14SB')
ff.promote_parameter('ligand_forcefield', promoted_name='ligand_ff', default='OpenFF_2.0.0')

# Production run
prod = ParallelMDNptCube("Production", title="Production")
prod.promote_parameter('time', promoted_name='prod_ns',
                       default=2,
                       description='Length of MD run in nanoseconds')
prod.promote_parameter('trajectory_interval', promoted_name='prod_trajectory_interval',
                       default=0.004,
                       description='Trajectory saving interval in ns')
prod.promote_parameter('hmr', promoted_name="HMR", title='Use Hydrogen Mass Repartitioning', default=True,
                       description='Give hydrogens more mass to speed up the MD')
prod.promote_parameter('md_engine', promoted_name='md_engine', default='OpenMM',
                       description='Select the MD Engine')
prod.set_parameters(reporter_interval=0.004)
prod.set_parameters(suffix='prod')

# Minimization
minimization = ParallelMDMinimizeCube('minComplex', title='Minimization')
minimization.modify_parameter(minimization.restraints, promoted=False, default="noh (ligand or protein)")
minimization.modify_parameter(minimization.restraintWt, promoted=False, default=5.0)
minimization.modify_parameter(minimization.steps, promoted=False, default=0)
minimization.set_parameters(center=True)
minimization.set_parameters(save_md_stage=True)
minimization.set_parameters(hmr=False)
minimization.promote_parameter("md_engine", promoted_name="md_engine")

# NVT simulation. Here the assembled system is warmed up to the final selected temperature
warmup = ParallelMDNvtCube('warmup', title='Warm Up')
warmup.set_parameters(time=0.01)
warmup.modify_parameter(warmup.restraints, promoted=False, default="noh (ligand or protein)")
warmup.modify_parameter(warmup.restraintWt, promoted=False, default=2.0)
warmup.set_parameters(trajectory_interval=0.0)
warmup.set_parameters(reporter_interval=0.001)
warmup.set_parameters(suffix='warmup')
warmup.set_parameters(hmr=False)
warmup.set_parameters(save_md_stage=True)
warmup.promote_parameter("md_engine", promoted_name="md_engine")

# The system is equilibrated at the right pressure and temperature in several stages
# The main difference between the stages is related to the restraint force used
# to keep the ligand and protein in their starting positions. A relatively strong force
# is applied in the first stage while a relatively small one is applied in the latter

# NPT Equilibration stage 1
equil1 = ParallelMDNptCube('equil1', title='Equilibration I')
equil1.set_parameters(time=0.01)
equil1.promote_parameter("hmr", promoted_name="HMR", default=True)
equil1.modify_parameter(equil1.restraints, promoted=False, default="noh (ligand or protein)")
equil1.modify_parameter(equil1.restraintWt, promoted=False, default=1.0)
equil1.set_parameters(trajectory_interval=0.0)
equil1.set_parameters(reporter_interval=0.001)
equil1.set_parameters(suffix='equil1')
equil1.promote_parameter("md_engine", promoted_name="md_engine")

# NPT Equilibration stage 2
equil2 = ParallelMDNptCube('equil2', title='Equilibration II')
equil2.set_parameters(time=0.02)
equil2.promote_parameter("hmr", promoted_name="HMR", default=True)
equil2.modify_parameter(equil2.restraints, promoted=False, default="noh (ligand or protein)")
equil2.modify_parameter(equil2.restraintWt, promoted=False, default=0.5)
equil2.set_parameters(trajectory_interval=0.0)
equil2.set_parameters(reporter_interval=0.001)
equil2.set_parameters(suffix='equil2')
equil2.promote_parameter("md_engine", promoted_name="md_engine")

# NPT Equilibration stage 3
equil3 = ParallelMDNptCube('equil3', title='Equilibration III')
equil3.modify_parameter(equil3.time, promoted=False, default=0.1)
equil3.promote_parameter("hmr", promoted_name="HMR", default=True)
equil3.modify_parameter(equil3.restraints, promoted=False, default="noh (ligand or protein)")
equil3.modify_parameter(equil3.restraintWt, promoted=False, default=0.2)
equil3.set_parameters(trajectory_interval=0.0)
equil3.set_parameters(reporter_interval=0.002)
equil3.set_parameters(suffix='equil3')
equil3.promote_parameter("md_engine", promoted_name="md_engine")

# NPT Equilibration stage 4
equil4 = ParallelMDNptCube('equil4', title='Equilibration IV')
equil4.modify_parameter(equil4.time, promoted=False, default=0.1)
equil4.promote_parameter("hmr", promoted_name="HMR", default=True)
equil4.modify_parameter(equil4.restraints, promoted=False, default="ca_protein or (noh ligand)")
equil4.modify_parameter(equil4.restraintWt, promoted=False, default=0.1)
equil4.set_parameters(trajectory_interval=0.0)
equil4.set_parameters(reporter_interval=0.002)
equil4.set_parameters(suffix='equil4')
equil4.promote_parameter("md_engine", promoted_name="md_engine")

md_group = ParallelCubeGroup(cubes=[minimization, warmup, equil1, equil2, equil3, equil4, prod])
job.add_group(md_group)


job.add_cubes(iligs, ligset, chargelig, ligid, iprot, mdcomp, complx, solvate,
              coll_open, ff, minimization, warmup, equil1, equil2, equil3, equil4, prod,
              coll_close, check_rec, exceptions, ofs, fail)

# Success Connections
iligs.success.connect(ligset.intake)
ligset.success.connect(chargelig.intake)
chargelig.success.connect(ligid.intake)
ligid.success.connect(complx.intake)
iprot.success.connect(mdcomp.intake)
mdcomp.success.connect(complx.protein_port)
complx.success.connect(solvate.intake)
solvate.success.connect(coll_open.intake)
coll_open.success.connect(ff.intake)

ff.success.connect(minimization.intake)
minimization.success.connect(warmup.intake)
warmup.success.connect(equil1.intake)
equil1.success.connect(equil2.intake)
equil2.success.connect(equil3.intake)
equil3.success.connect(equil4.intake)
equil4.success.connect(prod.intake)

prod.success.connect(coll_close.intake)
coll_close.success.connect(check_rec.intake)
check_rec.success.connect(ofs.intake)
exceptions.success.connect(fail.intake)

# Fail Connections
ligset.failure.connect(check_rec.fail_in)
chargelig.failure.connect(check_rec.fail_in)
ligid.failure.connect(check_rec.fail_in)
mdcomp.failure.connect(check_rec.fail_in)
complx.failure.connect(check_rec.fail_in)

ff.failure.connect(check_rec.fail_in)
minimization.failure.connect(check_rec.fail_in)
warmup.failure.connect(check_rec.fail_in)
equil1.failure.connect(check_rec.fail_in)
equil2.failure.connect(check_rec.fail_in)
equil3.failure.connect(check_rec.fail_in)
equil4.failure.connect(check_rec.fail_in)
prod.failure.connect(check_rec.fail_in)
coll_close.failure.connect(check_rec.fail_in)
check_rec.failure.connect(exceptions.intake)


if __name__ == "__main__":
    job.run()

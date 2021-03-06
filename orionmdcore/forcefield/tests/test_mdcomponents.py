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
    from datarecord import read_records
except ImportError:
    from orionmdcore import __installation__error__
    raise ImportError(__installation__error__)

from orionmdcore.standards import Fields

from orionmdcore.forcefield import MDComponents

import orionmdcore

import unittest

import os

PACKAGE_DIR = os.path.dirname(os.path.dirname(orionmdcore.__file__))
FILE_DIR = os.path.join(PACKAGE_DIR, "tests", "data")


class MDComponentsTest(unittest.TestCase):
    def test_Gaff_sdf(self):

        mol = oechem.OEMol()

        fn = os.path.join(FILE_DIR, "ex_mol.oeb")

        with oechem.oemolistream(fn) as ifs:
            oechem.OEReadMolecule(ifs, mol)

        md_components = MDComponents(from_molecule=mol)

        print(md_components)

        md_components.parametrize_components(
            protein_ff="Amber14SB", ligand_ff="Gaff_1.81", other_ff="Gaff_1.81"
        )

    def test_Gaff2(self):

        fn = os.path.join(FILE_DIR, "6puq_solvated.oedb")

        ifs = oechem.oeifstream(fn)

        records = []

        for rec in read_records(ifs):
            records.append(rec)
        ifs.close()

        record = records[0]

        md_components = record.get_value(Fields.md_components)

        print(md_components)

        md_components.parametrize_components(
            protein_ff="Amber14SB", ligand_ff="Gaff_2.11", other_ff="Gaff_2.11"
        )

    def test_Smirnoff99Frosst(self):

        fn = os.path.join(FILE_DIR, "6puq_solvated.oedb")

        ifs = oechem.oeifstream(fn)

        records = []

        for rec in read_records(ifs):
            records.append(rec)
        ifs.close()

        record = records[0]

        md_components = record.get_value(Fields.md_components)

        print(md_components)

        md_components.parametrize_components(
            protein_ff="Amber14SB",
            ligand_ff="Smirnoff99Frosst",
            other_ff="Smirnoff99Frosst",
        )

    def test_OpenFF1_1(self):

        fn = os.path.join(FILE_DIR, "6puq_solvated.oedb")

        ifs = oechem.oeifstream(fn)

        records = []

        for rec in read_records(ifs):
            records.append(rec)
        ifs.close()

        record = records[0]

        md_components = record.get_value(Fields.md_components)

        print(md_components)

        md_components.parametrize_components(
            protein_ff="Amber14SB", ligand_ff="OpenFF_1.1.1", other_ff="OpenFF_1.1.1"
        )

    def test_OpenFF1_2_1(self):
        fn = os.path.join(FILE_DIR, "6puq_solvated.oedb")

        ifs = oechem.oeifstream(fn)

        records = []

        for rec in read_records(ifs):
            records.append(rec)
        ifs.close()

        record = records[0]

        md_components = record.get_value(Fields.md_components)

        print(md_components)

        md_components.parametrize_components(
            protein_ff="Amber14SB", ligand_ff="OpenFF_1.2.1", other_ff="OpenFF_1.2.1"
        )

    def test_OpenFF1_3_1(self):
        fn = os.path.join(FILE_DIR, "6puq_solvated.oedb")

        ifs = oechem.oeifstream(fn)

        records = []

        for rec in read_records(ifs):
            records.append(rec)
        ifs.close()

        record = records[0]

        md_components = record.get_value(Fields.md_components)

        print(md_components)

        md_components.parametrize_components(
            protein_ff="Amber14SB", ligand_ff="OpenFF_1.3.1", other_ff="OpenFF_1.3.1"
        )

    def test_protein_force_field_amber_99sbildn_Gaff2(self):

        fn = os.path.join(FILE_DIR, "6puq_solvated.oedb")

        ifs = oechem.oeifstream(fn)

        records = []

        for rec in read_records(ifs):
            records.append(rec)
        ifs.close()

        record = records[0]

        md_components = record.get_value(Fields.md_components)

        print(md_components)

        md_components.parametrize_components(
            protein_ff="Amber99SBildn", ligand_ff="Gaff_2.11", other_ff="Gaff_2.11"
        )

    def test_protein_force_field_amber_99sbildn_OFF_1_2_1(self):

        fn = os.path.join(FILE_DIR, "6puq_solvated.oedb")

        ifs = oechem.oeifstream(fn)

        records = []

        for rec in read_records(ifs):
            records.append(rec)
        ifs.close()

        record = records[0]

        md_components = record.get_value(Fields.md_components)

        print(md_components)

        md_components.parametrize_components(
            protein_ff="Amber99SBildn",
            ligand_ff="OpenFF_1.2.1",
            other_ff="OpenFF_1.2.1",
        )

    def test_protein_force_field_amber_fb15_OFF_1_2_1(self):

        fn = os.path.join(FILE_DIR, "6puq_solvated.oedb")

        ifs = oechem.oeifstream(fn)

        records = []

        for rec in read_records(ifs):
            records.append(rec)
        ifs.close()

        record = records[0]

        md_components = record.get_value(Fields.md_components)

        print(md_components)

        md_components.parametrize_components(
            protein_ff="AmberFB15", ligand_ff="OpenFF_1.2.1", other_ff="OpenFF_1.2.1"
        )

    def test_protein_force_field_amber_99sbildn_OFF_1_3_1(self):

        fn = os.path.join(FILE_DIR, "6puq_solvated.oedb")

        ifs = oechem.oeifstream(fn)

        records = []

        for rec in read_records(ifs):
            records.append(rec)
        ifs.close()

        record = records[0]

        md_components = record.get_value(Fields.md_components)

        print(md_components)

        md_components.parametrize_components(
            protein_ff="Amber99SBildn",
            ligand_ff="OpenFF_1.3.1",
            other_ff="OpenFF_1.3.1",
        )

    def test_protein_force_field_amber_SB14_OFF_1_3_1(self):

        fn = os.path.join(FILE_DIR, "6puq_solvated.oedb")

        ifs = oechem.oeifstream(fn)

        records = []

        for rec in read_records(ifs):
            records.append(rec)
        ifs.close()

        record = records[0]

        md_components = record.get_value(Fields.md_components)

        print(md_components)

        md_components.parametrize_components(
            protein_ff="Amber14SB", ligand_ff="OpenFF_1.3.1", other_ff="OpenFF_1.3.1"
        )

    def test_protein_force_field_amber_fb15_OFF_1_3_1(self):

        fn = os.path.join(FILE_DIR, "6puq_solvated.oedb")

        ifs = oechem.oeifstream(fn)

        records = []

        for rec in read_records(ifs):
            records.append(rec)
        ifs.close()

        record = records[0]

        md_components = record.get_value(Fields.md_components)

        print(md_components)

        md_components.parametrize_components(
            protein_ff="AmberFB15", ligand_ff="OpenFF_1.3.1", other_ff="OpenFF_1.3.1"
        )

    def test_protein_force_field_amber_SB14_OFF_2_0_0(self):

        fn = os.path.join(FILE_DIR, "6puq_solvated.oedb")

        ifs = oechem.oeifstream(fn)

        records = []

        for rec in read_records(ifs):
            records.append(rec)
        ifs.close()

        record = records[0]

        md_components = record.get_value(Fields.md_components)

        print(md_components)

        md_components.parametrize_components(
            protein_ff="Amber14SB", ligand_ff="OpenFF_2.0.0", other_ff="OpenFF_2.0.0"
        )

    def test_protein_force_field_amber_fb15_OFF_2_0_0(self):

        fn = os.path.join(FILE_DIR, "6puq_solvated.oedb")

        ifs = oechem.oeifstream(fn)

        records = []

        for rec in read_records(ifs):
            records.append(rec)
        ifs.close()

        record = records[0]

        md_components = record.get_value(Fields.md_components)

        print(md_components)

        md_components.parametrize_components(
            protein_ff="AmberFB15", ligand_ff="OpenFF_2.0.0", other_ff="OpenFF_2.0.0"
        )

    def test_protein_force_field_amber_SB14_OFF_2_0_0_nc(self):

        fn = os.path.join(FILE_DIR, "Thr_3I.oedb")

        ifs = oechem.oeifstream(fn)

        records = []

        for rec in read_records(ifs):
            records.append(rec)
        ifs.close()

        record = records[0]

        md_components = record.get_value(Fields.md_components)

        print(md_components)

        md_components.parametrize_components(protein_ff='Amber14SB',
                                             ligand_ff='OpenFF_2.0.0',
                                             other_ff='OpenFF_2.0.0')


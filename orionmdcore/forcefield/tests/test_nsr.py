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

    from oeommtools.utils import oemol_to_openmmTop
except ImportError:
    from orionmdcore import __installation__error__
    raise ImportError(__installation__error__)

from orionmdcore.forcefield import nsr_template_generator

from orionmdcore.forcefield.ff_library import Default


from simtk.openmm.app import ForceField

from io import StringIO

import unittest


class NonStandardResidueTests(unittest.TestCase):
    """
    Test Remove Water and Ions
    """

    def test_nsr_no_termini_amber14(self):
        protein_fn = "tests/data/1h1q_no_termini.oeb"

        protein = oechem.OEMol()

        with oechem.oemolistream(protein_fn) as ifs:
            oechem.OEReadMolecule(ifs, protein)

        topology, position = oemol_to_openmmTop(protein)
        forcefield = ForceField("amber14/protein.ff14SB.xml")
        openff = Default.ligandff.offxml

        listff_templates = nsr_template_generator(
            protein, topology, forcefield, openff=openff
        )

        for ffxml_template in listff_templates:
            forcefield.loadFile(StringIO(ffxml_template))

    def test_nsr_no_termini_amber99(self):
        protein_fn = "tests/data/1h1q_no_termini.oeb"

        protein = oechem.OEMol()

        with oechem.oemolistream(protein_fn) as ifs:
            oechem.OEReadMolecule(ifs, protein)

        topology, position = oemol_to_openmmTop(protein)
        forcefield = ForceField("amber99sbildn.xml")
        openff = Default.ligandff.offxml

        listff_templates = nsr_template_generator(
            protein, topology, forcefield, openff=openff
        )

        for ffxml_template in listff_templates:
            forcefield.loadFile(StringIO(ffxml_template))

    def test_nsr_C_terminus_amber14(self):
        protein_fn = "tests/data/1c1b_C_termini.oeb"

        protein = oechem.OEMol()

        with oechem.oemolistream(protein_fn) as ifs:
            oechem.OEReadMolecule(ifs, protein)

        topology, position = oemol_to_openmmTop(protein)
        forcefield = ForceField("amber14/protein.ff14SB.xml")
        openff = Default.ligandff.offxml

        listff_templates = nsr_template_generator(
            protein, topology, forcefield, openff=openff
        )

        for ffxml_template in listff_templates:
            forcefield.loadFile(StringIO(ffxml_template))

    def test_nsr_C_terminus_amber99(self):
        protein_fn = "tests/data/1c1b_C_termini.oeb"

        protein = oechem.OEMol()

        with oechem.oemolistream(protein_fn) as ifs:
            oechem.OEReadMolecule(ifs, protein)

        topology, position = oemol_to_openmmTop(protein)
        forcefield = ForceField("amber99sbildn.xml")
        openff = Default.ligandff.offxml

        listff_templates = nsr_template_generator(
            protein, topology, forcefield, openff=openff
        )

        for ffxml_template in listff_templates:
            forcefield.loadFile(StringIO(ffxml_template))

    def test_nsr_N_terminus_amber14(self):
        protein_fn = "tests/data/1ly7_N_termini.oeb"

        protein = oechem.OEMol()

        with oechem.oemolistream(protein_fn) as ifs:
            oechem.OEReadMolecule(ifs, protein)

        topology, position = oemol_to_openmmTop(protein)
        forcefield = ForceField("amber14/protein.ff14SB.xml")
        openff = Default.ligandff.offxml

        listff_templates = nsr_template_generator(
            protein, topology, forcefield, openff=openff
        )

        for ffxml_template in listff_templates:
            forcefield.loadFile(StringIO(ffxml_template))

    def test_nsr_N_terminus_amber99(self):
        protein_fn = "tests/data/1ly7_N_termini.oeb"

        protein = oechem.OEMol()

        with oechem.oemolistream(protein_fn) as ifs:
            oechem.OEReadMolecule(ifs, protein)

        topology, position = oemol_to_openmmTop(protein)
        forcefield = ForceField("amber99sbildn.xml")
        openff = Default.ligandff.offxml

        listff_templates = nsr_template_generator(
            protein, topology, forcefield, openff=openff
        )

        for ffxml_template in listff_templates:
            forcefield.loadFile(StringIO(ffxml_template))

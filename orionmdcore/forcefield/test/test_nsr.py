from openeye import oechem

from orionmdcore.forcefield import nsr_template_generator

from oeommtools.utils import oemol_to_openmmTop

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
        openff = "openff_unconstrained-1.3.0.offxml"

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
        openff = "openff_unconstrained-1.3.0.offxml"

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
        openff = "openff_unconstrained-1.3.0.offxml"

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
        openff = "openff_unconstrained-1.3.0.offxml"

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
        openff = "openff_unconstrained-1.3.0.offxml"

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
        openff = "openff_unconstrained-1.3.0.offxml"

        listff_templates = nsr_template_generator(
            protein, topology, forcefield, openff=openff
        )

        for ffxml_template in listff_templates:
            forcefield.loadFile(StringIO(ffxml_template))

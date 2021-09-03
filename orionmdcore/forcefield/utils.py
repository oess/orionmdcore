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
    from oeommtools import utils as oeommutils

    from openeye import oechem, oequacpac, oeomega
except ImportError:
    from orionmdcore import __installation__error__
    raise ImportError(__installation__error__)

import numpy as np

import tempfile

import parmed

import openmoltools

import os

import shutil

from simtk.openmm.app import AmberInpcrdFile, AmberPrmtopFile

from simtk.openmm import app

from openff.toolkit.typing.engines.smirnoff import ForceField

from openff.toolkit.topology import Topology, Molecule

from orionmdcore.forcefield.ff_library import ff_library

import mdtraj.utils

import subprocess

from importlib.machinery import PathFinder


# TODO TEMPORARY SOLUTION FOR OPENMOLTOOLS BUG
# https://github.com/choderalab/openmoltools/issues/299


def run_tleap(
    molecule_name,
    gaff_mol2_filename,
    frcmod_filename,
    prmtop_filename=None,
    inpcrd_filename=None,
    leaprc="leaprc.gaff",
):

    """Run AmberTools tleap to create simulation files for AMBER

    Parameters
    ----------
    molecule_name : str
        The name of the molecule
    gaff_mol2_filename : str
        GAFF format mol2 filename produced by antechamber
    frcmod_filename : str
        Amber frcmod file produced by prmchk
    prmtop_filename : str, optional, default=None
        Amber prmtop file produced by tleap, defaults to molecule_name
    inpcrd_filename : str, optional, default=None
        Amber inpcrd file produced by tleap, defaults to molecule_name
    leaprc : str, optional, default = 'leaprc.gaff'
        Optionally, specify alternate leaprc to use, such as `leaprc.gaff2`

    Returns
    -------
    prmtop_filename : str
        Amber prmtop file produced by tleap
    inpcrd_filename : str
        Amber inpcrd file produced by tleap
    """
    if prmtop_filename is None:
        prmtop_filename = "%s.prmtop" % molecule_name
    if inpcrd_filename is None:
        inpcrd_filename = "%s.inpcrd" % molecule_name

    # Get absolute paths for input/output
    gaff_mol2_filename = os.path.abspath(gaff_mol2_filename)
    frcmod_filename = os.path.abspath(frcmod_filename)
    prmtop_filename = os.path.abspath(prmtop_filename)
    inpcrd_filename = os.path.abspath(inpcrd_filename)

    # Work in a temporary directory, on hard coded filenames,
    # to avoid any issues AMBER may have with spaces and other special characters in filenames
    with mdtraj.utils.enter_temp_directory():
        shutil.copy(gaff_mol2_filename, "file.mol2")
        shutil.copy(frcmod_filename, "file.frcmod")

        tleap_input = (
            """
    source oldff/leaprc.ff99SB
    source %s
    LIG = loadmol2 file.mol2
    loadamberparams file.frcmod
    check LIG
    saveamberparm LIG out.prmtop out.inpcrd
    quit

"""
            % leaprc
        )

        file_handle = open("tleap_commands", "w")
        file_handle.writelines(tleap_input)
        file_handle.close()

        cmd = "tleap -f %s " % file_handle.name

        subprocess.getoutput(cmd)

        # Copy back target files
        shutil.copy("out.prmtop", prmtop_filename)
        shutil.copy("out.inpcrd", inpcrd_filename)

    return prmtop_filename, inpcrd_filename


def assignELF10charges(molecule, max_confs=800, strictStereo=True, opt=None):
    """
     This function computes atomic partial charges for an OEMol by
     using the ELF10 method

    Parameters:
    -----------
    molecule : OEMol object
        The molecule that needs to be charged
    max_confs : integer
        The max number of conformers used to calculate the atomic partial charges
    strictStereo : bool
        a flag used to check if atoms need to have assigned stereo chemistry or not

    Return:
    -------
    mol_copy : OEMol
        a copy of the original molecule with assigned atomic partial charges
    """

    def generate_conformers(molecule, max_confs=800, strictStereo=True, ewindow=15.0, rms_threshold=1.0,
                            strictTypes=True):
        """Generate conformations for the supplied molecule

        Parameters
        ----------
        molecule : OEMol
            Molecule for which to generate conformers
        max_confs : int, optional, default=800
            Max number of conformers to generate.  If None, use default OE Value.
        strictStereo : bool, optional, default=True
            If False, permits smiles strings with unspecified stereochemistry.
        strictTypes : bool, optional, default=True
            If True, requires that Omega have exact MMFF types for atoms in molecule; otherwise, allows the closest atom type of the same element to be used.

        Returns
        -------
        molcopy : OEMol
            A multi-conformer molecule with up to max_confs conformers.

        Notes
        -----
        Roughly follows
        http://docs.eyesopen.com/toolkits/cookbook/python/modeling/am1-bcc.html

        """

        molcopy = oechem.OEMol(molecule)
        omega = oeomega.OEOmega()

        # These parameters were chosen to match http://docs.eyesopen.com/toolkits/cookbook/python/modeling/am1-bcc.html
        omega.SetMaxConfs(max_confs)
        omega.SetIncludeInput(True)
        omega.SetCanonOrder(False)

        omega.SetSampleHydrogens(
            True)  # Word to the wise: skipping this step can lead to significantly different charges!
        omega.SetEnergyWindow(ewindow)
        omega.SetRMSThreshold(
            rms_threshold)  # Word to the wise: skipping this step can lead to significantly different charges!

        omega.SetStrictStereo(strictStereo)
        omega.SetStrictAtomTypes(strictTypes)

        omega.SetIncludeInput(False)  # don't include input
        if max_confs is not None:
            omega.SetMaxConfs(max_confs)

        status = omega(molcopy)  # generate conformation
        if not status:
            raise (RuntimeError("omega returned error code %d" % status))

        return molcopy

    mol_copy_charged = oechem.OEMol(molecule)

    # The passed molecule could have already conformers. If the conformer number
    # does not exceed the max_conf threshold then max_confs conformations will
    # be generated
    if not mol_copy_charged.GetMaxConfIdx() > 200:
        # Generate up to max_confs conformers
        mol_copy_charged = generate_conformers(mol_copy_charged, max_confs=max_confs, strictStereo=strictStereo)

    # Assign MMFF Atom types
    if not oechem.OEMMFFAtomTypes(mol_copy_charged):
        raise RuntimeError("MMFF atom type assignment returned errors")

    # Check for Carboxylic Acid patterns in the molecule
    smarts = '(O=)[C][O,S][H]'
    ss = oechem.OESubSearch(smarts)

    oechem.OEPrepareSearch(mol_copy_charged, ss)
    unique_match = True

    a_match_list = []
    for match in ss.Match(mol_copy_charged, unique_match):

        for ma in match.GetAtoms():
            a_match_list.append(ma.target)

    # Set the Carboxylic Acid torsion to zero for each generated conformers
    if a_match_list:

        if len(a_match_list) % 4 != 0:
            raise ValueError("The atom matching list must be multiple of 4")

        for i in range(0, len(a_match_list), 4):

            chunk = a_match_list[i:i + 4]

            for conf in mol_copy_charged.GetConfs():

                conf.SetTorsion(chunk[0],
                                chunk[1],
                                chunk[2],
                                chunk[3], 0.0)

    # Try to calculate the ELF10 charges for the molecule
    quacpac_status = oequacpac.OEAssignCharges(mol_copy_charged, oequacpac.OEAM1BCCELF10Charges())

    if not quacpac_status:
        print("WARNING: OEAM1BCCELF10 charge assignment failed downgrading "
              "to OEAM1BCC charge assignment for the molecule: {}".format(mol_copy_charged.GetTitle()), flush=True)

        quacpac_status = oequacpac.OEAssignCharges(mol_copy_charged, oequacpac.OEAM1BCCCharges())

        if not quacpac_status:
            print("WARNING: OEAM1BCC charge assignment failed downgrading "
                  "to MMFF94 charge assignment for the molecule: {}".format(mol_copy_charged.GetTitle()), flush=True)

        quacpac_status = oequacpac.OEAssignCharges(mol_copy_charged, oequacpac.OEMMFF94Charges())

    if not quacpac_status:
        raise RuntimeError("OEAssignCharges returned error code {}".format(quacpac_status))

    # Copy back the charges to a molecule copy
    # without any conformers information
    map_charges = {at.GetIdx(): at.GetPartialCharge() for at in mol_copy_charged.GetAtoms()}

    mol_copy = oechem.OEMol(molecule)

    for at in mol_copy.GetAtoms():
        at.SetPartialCharge(map_charges[at.GetIdx()])

    # Check Formal vs Partial charges
    mol_copy_formal_charge = 0
    for at in mol_copy.GetAtoms():
        mol_copy_formal_charge += at.GetFormalCharge()

    mol_copy_partial_charge = 0.0
    for at in mol_copy.GetAtoms():
        mol_copy_partial_charge += at.GetPartialCharge()

    if abs(mol_copy_formal_charge - mol_copy_partial_charge) > 0.01:
        raise ValueError("Molecule Formal charge and Molecule Partial charge mismatch: {} vs {}".format(
            mol_copy_formal_charge, mol_copy_partial_charge))

    return mol_copy


class ParamMolStructure(object):
    """
    Generates parametrized ParmEd structure of the molecule with a chosen force field

    Parameters
    ----------
    molecule : openeye.oechem.OEMol
        The openeye molecule to be parameterized
    forcefield : str
        String specifying the forcefield parameters to be used
    prefix_name : str
        String specifying the output prefix filename

    Returns
    ---------
    packedmol : openeye.oechem.OEMol
        Openeye molecule with the ParmEd Structure attached.
    """

    def __init__(
        self,
        molecule,
        forcefield,
        prefix_name="ligand",
        delete_out_files=True,
        recharge=False,
    ):
        if forcefield not in list(ff_library.ligandff.values()):
            raise RuntimeError(
                "The selected ligand force field is not "
                "supported {}. Available {}".format(
                    forcefield, list(ff_library.ligandff.keys())
                )
            )
        else:
            self.molecule = molecule
            self.forcefield = str(forcefield).strip()
            self.structure = None
            self.prefix_name = prefix_name
            self.delete_out_files = delete_out_files
            self.recharge = recharge

    def checkCharges(self, molecule):
        # Check that molecule is charged.
        is_charged = False
        for atom in molecule.GetAtoms():
            if atom.GetPartialCharge() != 0.0:
                is_charged = True

        if is_charged and self.recharge:
            is_charged = False

        return is_charged

    def getSmirnoffStructure(self, molecule=None):

        if not molecule:
            molecule = self.molecule

        if not self.checkCharges(molecule):
            print("WARNING: Missing Charges, assigning elf10 charges to molecule")
            molecule = assignELF10charges(molecule)

        if self.forcefield == ff_library.ligandff["Smirnoff99Frosst"]:

            fffn = os.path.join(PathFinder().find_spec("openff").submodule_search_locations._path[0], "toolkit/data/test_forcefields/" + self.forcefield)

            # fffn = resource_filename(
            #     "openff-toolkit",
            #     os.path.join("toolkit", "data", "test_forcefields/" + self.forcefield),
            # )

            if not os.path.exists(fffn):
                raise ValueError(
                    "Sorry! {} does not exist. If you just added it, you'll have to re-install".format(
                        fffn
                    )
                )

            with open(fffn) as ffxml:
                ff = ForceField(ffxml, allow_cosmetic_attributes=True)

        elif self.forcefield in [v for k, v in ff_library.ligandff.items() if k not in ["Smirnoff99Frosst",
                                                                                        "Gaff_1.81",
                                                                                        "Gaff_2.11"]]:

            ff = ForceField(self.forcefield, allow_cosmetic_attributes=True)

        else:
            raise ValueError("Force Field not Supported: {}".format(self.forcefield))

        # Create a copy of the molecule otherwise the the following OFF function will
        # change the bonds as side effect. This was checked by using the Isomeric Smiles
        # molecule representation
        copy_mol = oechem.OEMol(molecule)

        mol_off = Molecule.from_openeye(copy_mol, allow_undefined_stereo=True)

        topology = Topology.from_molecules([mol_off])

        omm_sys = ff.create_openmm_system(topology, charge_from_molecules=[mol_off])

        # omm_top = generateTopologyFromOEMol(molecule)
        # positions = mol_off.conformers[0]

        omm_top, positions = oeommutils.oemol_to_openmmTop(molecule)

        pmd_structure = parmed.openmm.load_topology(omm_top, omm_sys, xyz=positions)

        return pmd_structure

    def getGaffStructure(self, molecule=None, forcefield=None):
        if not molecule:
            molecule = self.molecule

        if not forcefield:
            forcefield = self.forcefield

        if not molecule:
            molecule = self.molecule

        if not self.checkCharges(molecule):
            print("WARNING: Missing Charges, assigning elf10 charges to molecule")
            molecule = assignELF10charges(molecule)

        # ofs = oechem.oemolostream("pippo.mol2")
        # mol2_format = oechem.OEFormat_MOL2
        # flavor = oechem.OEOFlavor_MOL2_GeneralFFFormat
        # ofs.SetFormat(mol2_format)
        # ofs.SetFlavor(mol2_format, flavor)
        # oechem.OEClearAromaticFlags(molecule)
        # oechem.OETriposAtomTypeNames(molecule)
        # oechem.OETriposAtomNames(molecule)
        # oechem.OETriposBondTypeNames(molecule)
        # oechem.OEWriteMol2File(ofs, molecule)
        # ofs.close()

        mol_sdf_file = tempfile.NamedTemporaryFile(suffix=".sdf")
        mol_sdf_file = mol_sdf_file.name

        with oechem.oemolostream(mol_sdf_file) as ofs:
            oechem.OEWriteConstMolecule(ofs, molecule)

        prefix = self.prefix_name + "_" + os.path.basename(mol_sdf_file).split(".")[0]

        total_formal_charge = 0
        for at in molecule.GetAtoms():
            total_formal_charge += at.GetFormalCharge()

        gaff_mol2_filename, frcmod_filename = openmoltools.amber.run_antechamber(
            prefix,
            mol_sdf_file,
            input_format="sdf",
            gaff_version=forcefield.lower(),
            charge_method="gas",
            net_charge=total_formal_charge,
        )
        prmtop, inpcrd = run_tleap(
            prefix,
            gaff_mol2_filename,
            frcmod_filename,
            leaprc="leaprc.{}".format(forcefield.lower()),
        )

        # shutil.copy(prmtop, "./debug")

        # TODO Load via ParmEd: This is causing Problems
        #  Merging two structures (OpenMM PMD structure and
        #  Amber PMD Structure): The NB exception list is messed up
        # molecule_structure = parmed.amber.AmberParm(prmtop, inpcrd)

        # TODO MODIFIED BY GAC
        omm_prmtop = AmberPrmtopFile(prmtop)
        omm_inpcrd = AmberInpcrdFile(inpcrd)

        omm_system = omm_prmtop.createSystem(nonbondedMethod=app.NoCutoff)

        molecule_structure = parmed.openmm.load_topology(
            omm_prmtop.topology, omm_system, xyz=omm_inpcrd.positions
        )

        for atmol in molecule.GetAtoms():
            res = oechem.OEAtomGetResidue(atmol)
            res_name = res.GetName()
            break

        for at in molecule_structure.atoms:
            at.residue.name = res_name

        if molecule.NumAtoms() != len(molecule_structure.atoms):
            raise ValueError("OE and Parmed Molecule size mismatch")

        for oe_at, pmd_at in zip(molecule.GetAtoms(), molecule_structure.atoms):
            if oe_at.GetAtomicNum() != pmd_at.atomic_number:
                raise ValueError(
                    "Atomic number mismatch between Parmed and the OpenEye topologies: {} - {}".format(
                        oe_at.GetAtomicNum(), pmd_at.atomic_number
                    )
                )

            pmd_at.charge = oe_at.GetPartialCharge()

        if self.delete_out_files:
            os.remove(gaff_mol2_filename)
            os.remove(frcmod_filename)
            os.remove(prmtop)
            os.remove(inpcrd)

        return molecule_structure

    def parameterize(self):

        ligandff = [v for k, v in ff_library.ligandff.items() if k not in ["Gaff_1.81", "Gaff_2.11"]]

        if self.forcefield in ligandff:

            structure = self.getSmirnoffStructure()

        elif self.forcefield in ["GAFF", "GAFF2"]:
            structure = self.getGaffStructure()

        self.structure = structure
        return self.structure


def parametrize_component(component, component_ff, other_ff):

    component_copy = oechem.OEMol(component)

    # OpenMM topology and positions from OEMol
    topology, positions = oeommutils.oemol_to_openmmTop(component_copy)

    # Try to apply the selected FF on the component
    if isinstance(component_ff, list):
        forcefield = app.ForceField()
        for f in component_ff:
            forcefield.loadFile(f)
    else:
        forcefield = app.ForceField(component_ff)

    # List of the unrecognized component
    unmatched_res_list = forcefield.getUnmatchedResidues(topology)

    # Try to parametrize the whole flask
    if not unmatched_res_list:
        omm_components = forcefield.createSystem(
            topology, rigidWater=False, constraints=None
        )
        components_pmd = parmed.openmm.load_topology(
            topology, omm_components, xyz=positions
        )
        return components_pmd

    # Extract The different non bonded parts
    numparts, parts = oechem.OEDetermineComponents(component_copy)
    pred = oechem.OEPartPredAtom(parts)

    part_mols = dict()
    for i in range(1, numparts + 1):
        pred.SelectPart(i)
        partmol = oechem.OEMol()
        oechem.OESubsetMol(partmol, component_copy, pred)
        part_mols[i] = partmol

    part = -1
    part_numbers = []
    for at in component_copy.GetAtoms():
        if parts[at.GetIdx()] != part:
            part = parts[at.GetIdx()]
            part_numbers.append(part)
        # print("atom %d is in part %d" % (at.GetIdx(), parts[at.GetIdx()]))

    part_numbers_resolved = dict()
    for i in range(0, len(part_numbers)):
        if part_numbers[i] in part_numbers_resolved:
            continue
        else:
            part_numbers_resolved[part_numbers[i]] = part_numbers[i]
            if i == len(part_numbers) - 1:
                break
            for j in range(i + 1, len(part_numbers)):
                if part_numbers[j] in part_numbers_resolved:
                    continue
                else:
                    smiles_i = oechem.OECreateCanSmiString(part_mols[part_numbers[i]])
                    smiles_j = oechem.OECreateCanSmiString(part_mols[part_numbers[j]])

                    if smiles_i == smiles_j:
                        part_numbers_resolved[part_numbers[j]] = part_numbers[i]
                    else:
                        continue

    ordered_parts = [part_numbers_resolved[pn] for pn in part_numbers]

    # print(ordered_parts)

    unique_parts = set(ordered_parts)

    part_to_pmd = dict()
    for part_i in unique_parts:

        mol_part_i = part_mols[part_i]
        omm_top_part_i, omm_pos_part_i = oeommutils.oemol_to_openmmTop(mol_part_i)

        if not forcefield.getUnmatchedResidues(omm_top_part_i):

            omm_sys_part_i = forcefield.createSystem(
                omm_top_part_i, rigidWater=False, constraints=None
            )

            pmd_part_i = parmed.openmm.load_topology(
                omm_top_part_i, omm_sys_part_i, xyz=omm_pos_part_i
            )
            # print("Parametrize full part: {}".format(part_i))
        else:

            mol_part_i_clean = oeommutils.sanitizeOEMolecule(mol_part_i)
            pmd_part_i = ParamMolStructure(
                mol_part_i_clean,
                other_ff,
                prefix_name="Part" + "_" + str(part_i),
                recharge=True,
            ).parameterize()
            # print("Parametrize part: {}".format(part_i))
        part_to_pmd[part_i] = pmd_part_i

    # Component Parmed Structure
    component_pmd = parmed.Structure()
    for part_i in ordered_parts:
        component_pmd += part_to_pmd[part_i]

    # Checking
    for parm_at, oe_at in zip(component_pmd.atoms, component_copy.GetAtoms()):

        if parm_at.atomic_number != oe_at.GetAtomicNum():
            raise ValueError(
                "Atomic number mismatch between the Parmed and the OpenEye topologies: {} - {}".format(
                    parm_at.atomic_number, oe_at.GetAtomicNum()
                )
            )

    # Set the positions
    oe_comp_coord_dic = component_copy.GetCoords()
    comp_coords = np.ndarray(shape=(component_copy.NumAtoms(), 3))
    for at_idx in oe_comp_coord_dic:
        comp_coords[at_idx] = oe_comp_coord_dic[at_idx]

    component_pmd.coordinates = comp_coords

    return component_pmd


def clean_tags(molecule):
    """
    This function remove tags that could cause problems along the MD Analysis stage.
    In particular Hint interactions, Style and PDB data are removed.

    Parameters:
    -----------
    molecule: OEMol molecule
        The molecule to clean

    Return:
    -------
    molecule: OEMol molecule
        The cleaned molecule
    """

    oechem.OEDeleteInteractionsHintSerializationData(molecule)
    oechem.OEDeleteInteractionsHintSerializationIds(molecule)
    oechem.OEClearStyle(molecule)
    oechem.OEClearPDBData(molecule)

    return molecule

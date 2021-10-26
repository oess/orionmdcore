* Purpose:

  * This Floe performs MD simulations given a prepared protein and a set of posed
    and prepared ligands.
* Method Recommendations/Requirements:

  * The ligands need to have reasonable 3D coordinates, all atoms, and correct
    chemistry (in particular bond orders and formal charges).
  * Each ligand can have multiple conformers but each conformer will be run
    separately as a different ligand.
  * The starting poses should not have very high gradients, in particular
    no bad clashes with the protein.
  * The protein needs to be prepared to MD standards: protein chains must be
    capped, all atoms in protein residues (including hydrogens) must be present,
    and missing protein loops resolved or capped.
  * Crystallographic internal waters should be retained where possible.
* Limitations

  * Currently this floe cannot handle covalent bonds between different components
    such as ligand, protein, and cofactors.
  * Glycosylation on proteins is truncated and the amino acid is capped with H.
* Expertise Level:

  * Regular/Intermediate/Advanced
* Compute Resource:

  * Depends on simulation length; Minimal resources for default 2 ns.
* Keywords:

  * MD, MDPrep
* Related Floes:

  * Short Trajectory MD with Analysis [MDPrep] [MD]
  * Ligand Bound and Unbound Equilibration for NES [MDPrep] [MD]

Given the inputs of the protein and posed ligands,
the complex is formed with each ligand/conformer
separately, and the complex is solvated and parametrized according to the
selected force fields. A minimization stage is performed on the system followed
by a warm up (NVT ensemble) and several equilibration stages (NPT ensemble). In the
minimization, warm up, and equilibration stages, positional harmonic restraints are
applied on the ligand and protein. At the end of the equilibration stages a 
production run (by default only 2 ns) is performed on the unrestrained system.
The trajectory and final state are written out in the results record;
no analysis is performed.

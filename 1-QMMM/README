## Workflow Overview

Starting from a PDB bound to ATP reactants:

1. Homology modeling (SWISS-MODEL): incorporate catalytic waters before adding the solvent box and RNA, generating `paratleap.pdb`.
2. Use `tleap` to solvate and add missing atoms, producing `x.prmtop` and `x.pdb` (reactants).
3. Perform minimization, heating, and equilibration.
4. Move atoms to generate `P.pdb` (optimized with ORCA prior to NEB). The topology must remain identical: same `x.prmtop`.

### NEB runs with ORCA (B3LYP 6‑31+G*)

1. Starting from `x.prmtop` and `R.pdb` (reactants), convert the topology for ORCA:
/path/to/orca/orca_mm -convff -AMBER x.prmtop
This generates `x.ORCAFF.prms`.
2. Run `NEB.inp` with optimization enabled (=True).  
On clusters, runs may be interrupted due to time limits; restart using the `restart` folder and the modified `NEB.inp`, copying necessary files from the interrupted run.
3. NEB has an internal convergence criterion and will terminate automatically (see ORCA manual).
4. `NEB.out` contains energy values; `NEB_MEP_trj.xyz` contains the converged band images.

### ASM runs with AMBER (DFTB3)

1. `NEB_reactant.pdb` and `NEB_product.pdb` (or images 1 and 15 from `NEB_MEP_trj.xyz`) are starting points for reactant and product coordinates.
2. Use `parmed` to replace H atoms in the QM region with D.
3. Minimize the structure.
4. Run a QM/MM dynamics equilibration and evaluate stability of reactants and products. These coordinates are used for string method input.
5. Choose collective variables and prepare the CVs file.
6. From the chosen CVs, calculate their values in `NEB_MEP_trj.xyz` to generate the guess file.

### References and Tutorials

- https://carlosramosg.com/category/tutorials  
- https://github.com/kzinovjev/string-amber/blob/master/tutorial.md  
- https://pubs.acs.org/doi/full/10.1021/acs.jpca.7b10842

# Ribose-binding protein: Cartesian vs torsion normal modes

A runnable, target-guided **backbone structural fitting benchmark**, with bundled source structures, numerical checks, measured outputs and an offline interactive viewer.

## Start here

1. Open `results/viewer.html` in a modern browser. No installation or network is needed to view the saved runs. Drag to rotate, scroll to zoom, select the mode count and trial, and play the amplitude interpolation.
2. Read `results/RESULTS.md` for measured outcomes and limitations.
3. To reproduce in Python 3.10 or later:

```bash
python -m venv .venv
```

Activate the environment on Windows PowerShell:

```powershell
.venv\Scripts\Activate.ps1
```

Or on macOS/Linux:

```bash
source .venv/bin/activate
```

Then:

```bash
python -m pip install -r requirements.txt
python check_model.py
python benchmark.py
python report.py
```

A smaller run:

```bash
python benchmark.py --modes 2 10 --seeds 0 --maxiter 150
python report.py
```

Running again overwrites matching result files and regenerates the manifest, summary and viewer for that invocation. Other run JSON files may remain; the report reads only runs listed in the current manifest. `report.py` requires matplotlib; the benchmark and viewer builder require only NumPy and SciPy. Record exact installed versions in `results/manifest.json`; requirements allow compatible versions, not bitwise-identical replay.

## Scope and inputs

- Start: 1BA2 chain A, D67R mutant, open, ligand-free.
- Reference: 2DRI chain A, wild type, closed, ribose-bound.
- Source PDB files are bundled unchanged; SHA-256 hashes are in the manifest.
- Modeled atoms: N, CA and C from each of 271 consecutive residues, 813 atoms total.
- At residue 67, the common backbone is matched despite ARG/ASP identity differences.
- Waters, ligand, side chains, oxygen atoms and hydrogens are excluded.
- The sequence/ligand differences remain scientific limitations even when only common backbone atoms are used.
- This is not an exact reproduction of Bray et al.'s all-atom modes or their optimizer. It is an explicitly specified first implementation using their structural pair and fractional RMSD idea.

## Model and equations

Distances are in Å and angles in radians. The spring strengths define relative model stiffness, not calibrated molecular energies or physical time.

The same backbone elastic network is used for both methods:

$$
U(x)=\frac12\sum_{(i,j)\in C}\gamma_{ij}(\|x_i-x_j\|-d^A_{ij})^2.
$$

Pairs within 8 Å in the reference A form the fixed network. Spring strength is 10 for consecutive backbone atoms, 5 for separation by two backbone indices, and 1 otherwise. These are declared prototype choices, not fitted force-field parameters. Sensitivity to cutoff and weights should be studied before generalizing results. A is the network's zero-energy minimum by construction.

Cartesian modes solve:

$$
H_x e_k=\lambda_k e_k,\qquad H_x=\nabla^2 U(x_A).
$$

The six rigid-body modes are excluded. Candidates are:

$$
x_C(a)=x_A+E_K a.
$$

Torsions are actual backbone phi and psi rotations, implemented as sequential Rodrigues rotations of downstream backbone atoms. All peptide omega angles and proline phi angles stay fixed, as do backbone bond lengths and bond angles. The first phi is excluded as a pure overall rotation and the last psi is absent. There are 531 selected torsions for this backbone. Rings/side chains are not modeled.

At the stationary reference:

$$
J_0=\frac{\partial f}{\partial q},\quad J=(I-QQ^T)J_0,
\quad H_q=J^T H_x J,\quad G=J^T J.
$$

Q is an orthonormal basis for overall translation and rotation of the modeled backbone. Solve:

$$
H_q v_k=\lambda_k Gv_k,\qquad x_T(b)=f(q_A+V_Kb).
$$

Metric eigenvalues below 1e-10 of the maximum are removed before whitening; retained rank is reported. This is an equal-atom-weight geometric stiffness comparison, not a mass-weighted vibrational-frequency calculation. Both mode bases are rescaled to have a unit infinitesimal Cα RMS displacement per coefficient; Cartesian vectors are not required to remain orthonormal after this scaling.

The target-guided objective for both methods is:

$$
\min_a\operatorname{RMSD}_{CA}(x_C(a),x_B)^2,
\qquad
\min_b\operatorname{RMSD}_{CA}(x_T(b),x_B)^2.
$$

Every objective evaluation optimally aligns the structures by a proper Kabsch rotation on Cα atoms. Both branches use L-BFGS-B, analytic gradients, mode-amplitude bounds [-8, 8], maximum 150 iterations, ftol=1e-10 and gtol=1e-6. Amplitudes have been normalized by tangent Cα displacement; equal coefficients do not imply equal finite physical motion in nonlinear torsion space. Seed 0 starts at zero; seeds 1 and 2 start with Gaussian coefficients of standard deviation 0.25. These bounds and initializations are prototype settings. Target B guides the fit, so this is not blind prediction.

$$
f_{RMS}=\frac{\operatorname{RMSD}_{CA}(X,B)}{\operatorname{RMSD}_{CA}(A,B)}.
$$

The mode counts are 1, 2, 5 and 10. Modes are computed from the fixed starting structure, not relinearized along the motion. There is no physical energy ranking, dynamics integration, free-energy estimate or geometry relaxation.

## Fairness and timing

Each method/mode-count/seed runs in a fresh, sequential subprocess with BLAS/OpenMP threads set to 1. Total measured compute time includes structure preparation, dense Hessian construction, eigensolution (including the torsion metric and Jacobian), optimization, final reconstruction, and final geometry checks. Display-frame generation and file writing are excluded and the former is timed separately. Setup is recomputed for every run, rather than silently amortized through a cache.

Peak RSS includes the Python interpreter, loaded libraries and frame generation as well as science computation; on platforms without resource it is null. It does not isolate the matrices' incremental memory. All Cartesian runs precede torsion runs; the pilot is not a randomized-order performance study. Three seeds are initial repeatability evidence, not confidence intervals for general algorithm performance.

The accuracy-time viewer starts from the baseline and shows best evaluated RMSD versus elapsed computation. Optimizer trial points are included. Geometry is screened only for final candidates. The thresholds CSV therefore measures RMSD-only attainment, not attainment of validated molecular geometry. Runs not reaching a threshold are explicitly marked, never assigned made-up times.

Cartesian and torsion branches differ in constraints as well as coordinate representation. Their comparison is between practical implemented modeling pipelines; isolating coordinate-system overhead would require matching the feasible physical space and constraints.

## Validation

`check_model.py` checks:

- matched atom counts, residue order and the known sequence mismatch;
- rigid-body invariance of aligned RMSD;
- exact zero-angle reconstruction;
- preservation of bond lengths and angles under torsion moves;
- analytic torsion gradients against centered finite differences;
- Hessian symmetry and translational nullspace;
- Hessian curvature against finite differences of the independently evaluated network energy.

The Cartesian mode solver additionally requires six zero modes. Input chain breaks and incomplete backbone records fail explicitly.

Geometry outputs report maximum deviations from the starting backbone and the count of nonlocal backbone atom pairs below 2 Å (index separation at least four). This is a simple diagnostic, not MolProbity or full all-atom validation. Distorted Cartesian fits should be treated as geometric fits requiring further refinement, not usable molecular conformers.

## Project files

- `model.py`: structure parsing, alignment, network, torsion kinematics, modes and geometry.
- `benchmark.py`: isolated computation, fitting, histories and metadata.
- `check_model.py`: independent numerical checks.
- `viewer.py`, `viewer_template.html`: standalone result viewer generator.
- `report.py`: plots, aggregate measurements and results narrative.
- `data/`: original PDB files.
- `results/`: complete measured example run.

## References

- [1BA2, RCSB PDB](https://www.rcsb.org/structure/1BA2), DOI: 10.2210/pdb1BA2/pdb.
- [2DRI, RCSB PDB](https://www.rcsb.org/structure/2DRI), DOI: 10.2210/pdb2DRI/pdb.
- Bray, Weiss and Levitt (2011), *Optimized Torsion-Angle Normal Modes Reproduce Conformational Changes More Accurately Than Cartesian Modes*. DOI: [10.1016/j.bpj.2011.10.054](https://doi.org/10.1016/j.bpj.2011.10.054).
- Hinsen, *Normal mode theory and harmonic potential approximations*, supplied reference: single-well harmonic framework and its limits.
- [MDAnalysis BAT documentation, version 2.6.0](https://docs.mdanalysis.org/2.6.0/documentation_pages/analysis/bat.html): internal-coordinate reference. This implementation uses explicit backbone rotations and does not depend on MDAnalysis.

The supplied papers are referenced but not redistributed in this project.

## Thirteen protein benchmark

The `benchmark_13.py` runner covers all thirteen structure pairs in Bray, Weiss and Levitt (2011), with ten modes and one initialization seed in **both** directions. It adds an explicit relative six-degree-of-freedom pose for every additional matched chain. See [`results_13/REPORT.md`](results_13/REPORT.md) and open `results_13/index.html` for the interactive per-protein viewers.

```bash
python benchmark_13.py --both-directions --k 10 --seed 0
```

This command runs every requested method and direction in a fresh, single-threaded process and writes a record for successes, failures or stopped optimizers. At ten modes and seed zero it generates the report and viewers. `--pairs 0 1` runs just the first two rows; pair indices begin at zero. Source PDBs and exact SHA-256 hashes are bundled. The original ribose run remains under `results/`.

**Matching rules.** The new runner matches complete N–CA–C residues by author residue number when at least 75% of the smaller protein matches and the overlapping residue identity is at least 80%. If residue numbering differs, it takes exact sequence blocks covering at least 75%. It reports the coverage, matching method and substitutions. Inputs with too few matches, incomplete specified chains, fewer than 15 shared residues, or more than 3500 modeled backbone atoms fail with an explanation. Ambiguous homologs or reordered chains require manual residue alignment before benchmarking.

**Internal backbone breaks.** A C–N gap exceeding 1.8 Å between matched residues creates separate segments. Torsions cannot rotate across missing backbone; the Cartesian network may have disconnected or weakly connected fragments. Five of the thirteen named pairs have one or more internal breaks in this implementation (calmodulin, diphtheria toxin, scallop myosin II, T7 RNA polymerase, and nitrogen regulatory protein C). This makes their comparison approximate. Their rows are included and explicitly flagged; do not cite them as fully modeled chains. The other eight pairs have no internal matched-chain breaks. Missing ligands, domains, contacts, and side chains may still affect the eight continuous pairs.

**Complexes.** Chain IDs in `pairs.json` give the correspondence. The selected extra chains each have three additional translations and three rotations, letting a complex change relative chain pose. A single rigid-body alignment is used to score the whole complex. This is a backbone-only structural approximation; it does not simulate binding or unbinding thermodynamics.

## Benchmark a new structure pair

Run the code from the project directory, with two local PDB files for conformations of the *same or closely related protein*:

```bash
python benchmark_13.py \
  --start /path/to/open.pdb --target /path/to/closed.pdb \
  --start-chains A --target-chains A \
  --name 'My protein' --output results_custom
```

Use corresponding chain letters such as `--start-chains AB --target-chains XY` for a two-chain complex. At ten modes, seed zero, `results_custom/index.html` and `results_custom/protein_01.html` provide an offline viewer, and `REPORT.md` and `summary.csv` give the measurements. `--both-directions` additionally runs target→start. Only local **PDB format** files are supported; FASTA sequences alone, mmCIF, arbitrary unrelated proteins, and all-atom physics are outside this implementation.

For several pairs, provide a JSON list in the same schema as `pairs.json` and run:

```bash
python benchmark_13.py --config my_pairs.json --output results_custom --both-directions
```

Every pair must contain `name`, `a`, `b`, `chains_a`, and `chains_b`. Set `a` and `b` to local PDB paths or IDs of PDB files already bundled in `data/`. Run with a new output directory for each dataset. The cutoff (8 Å), contact weights (10/5/1), mode count (10), and amplitude bounds (±8) were fixed for this first cohort. The model does not choose appropriate parameters automatically for an unfamiliar protein. Measure sensitivity and input quality before applying any biological interpretation.

## Interpretation of this cohort

The main score is target-guided aligned Cα fractional RMSD. A lower number means a better fit to the *known target structure*. Method rankings include differences in allowed motions, chain handling and missing-backbone treatment. Timing is specific to this machine and includes mode calculation; it is not a hardware-independent computational complexity result. The per-protein viewers animate mode amplitudes and **do not display molecular dynamics time**.

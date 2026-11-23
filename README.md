# Water benchmark: how basis set and method change computed properties

A small, reproducible PySCF study of one molecule. It asks two questions:

1. How do the **basis set** and its type (polarisation, diffuse functions) change the Hartree-Fock energy and dipole moment of water?
2. How do **Hartree-Fock (HF), PBE and B3LYP** compare with experiment for water's geometry and dipole moment, once each method has relaxed the geometry itself, and does that comparison depend on the basis set?

## Method

| Part | What is computed | Level of theory |
|---|---|---|
| A. Basis-set study | Energy and dipole at a fixed, near-experimental geometry (O-H 0.958 A, 104.4 deg) | HF with STO-3G, 6-31G, 6-31G*, cc-pVDZ, cc-pVTZ, aug-cc-pVDZ, aug-cc-pVTZ |
| B. Method study | Optimised O-H bond length, H-O-H angle and dipole | HF, PBE, B3LYP, each with cc-pVTZ and aug-cc-pVTZ |

- Software: [PySCF](https://pyscf.org) for the electronic structure, [geomeTRIC](https://github.com/leeping/geomeTRIC) for geometry optimisation.
- Pople basis sets use Cartesian d functions (their standard definition, so 6-31G* has 19 functions); Dunning sets use spherical functions.
- Experimental references (gas-phase water, NIST CCCBDB): r(OH) = 0.9572 A, angle = 104.52 deg, dipole = 1.855 D.

## Results

### A. Basis-set study (HF, fixed geometry)

Energy drops are only given within a basis family, where each set is a systematic improvement on the previous one.

| family | basis | basis_functions | energy_Eh | energy_drop_kcal_mol | dipole_D | dipole_error_D |
|---|---|---|---|---|---|---|
| Pople | sto-3g | 7 | -74.963063 |  | 1.726 | -0.129 |
| Pople | 6-31g | 13 | -75.983948 | 640.6 | 2.632 | +0.777 |
| Pople | 6-31g* | 19 | -76.010496 | 16.7 | 2.227 | +0.372 |
| Dunning | cc-pvdz | 24 | -76.026766 |  | 2.058 | +0.203 |
| Dunning | cc-pvtz | 58 | -76.057114 | 19.0 | 2.027 | +0.172 |
| aug-Dunning | aug-cc-pvdz | 41 | -76.041382 |  | 2.001 | +0.146 |
| aug-Dunning | aug-cc-pvtz | 92 | -76.060559 | 12.0 | 1.985 | +0.130 |

![Dipole versus basis set](results/dipole_vs_basis.png)

### B. Method study (geometry optimised for each method and basis)

| basis | method | r_OH_A | r_OH_error_A | angle_deg | angle_error_deg | dipole_D | dipole_error_D |
|---|---|---|---|---|---|---|---|
| cc-pvtz | HF | 0.9406 | -0.0166 | 106.00 | +1.48 | 1.988 | +0.133 |
| cc-pvtz | PBE | 0.9697 | +0.0125 | 103.53 | -0.99 | 1.885 | +0.030 |
| cc-pvtz | B3LYP | 0.9613 | +0.0041 | 104.54 | +0.02 | 1.919 | +0.064 |
| aug-cc-pvtz | HF | 0.9410 | -0.0162 | 106.32 | +1.80 | 1.939 | +0.084 |
| aug-cc-pvtz | PBE | 0.9703 | +0.0131 | 104.17 | -0.35 | 1.805 | -0.050 |
| aug-cc-pvtz | B3LYP | 0.9618 | +0.0046 | 105.10 | +0.58 | 1.847 | -0.008 |

![Error versus experiment](results/method_errors.png)

## Findings

- **The minimal basis set is not good enough.** Going from STO-3G to 6-31G lowers the energy by about 641 kcal/mol, far more than any chemical effect. STO-3G's dipole looks close to experiment only through fortuitous error cancellation, not because the basis is accurate.
- **Polarisation functions matter.** Adding d functions on oxygen (6-31G to 6-31G*) lowers the energy by about 17 kcal/mol and roughly halves the dipole error (+0.777 D to +0.372 D).
- **Diffuse functions matter more for the dipole than the zeta level.** aug-cc-pVDZ (41 functions) gives a smaller dipole error (+0.146 D) than cc-pVTZ (58 functions, +0.172 D). Going from cc-pVTZ to aug-cc-pVTZ lowers the energy by only about 2 kcal/mol but shifts the dipole by 0.04 D.
- **At the Hartree-Fock level the dipole converges to a method error.** It levels off near 2.0 D (1.985 D at aug-cc-pVTZ), about 0.13 D above experiment. The remaining error comes from the method, not the basis set. Relaxing the HF geometry (shorter O-H bond) lowers it further to 1.939 D.
- **Bond lengths are stable across the two basis sets; angles and dipoles are not.** HF gives bonds that are too short (-0.016 A). PBE gives bonds that are too long (+0.013 A) and B3LYP slightly long (+0.004 to +0.005 A). Bond lengths change by less than 0.001 A between cc-pVTZ and aug-cc-pVTZ, so B3LYP is closest for the bond in both.
- **The method ranking for the dipole and the angle flips with the basis set.**

  | property | cc-pVTZ ranking (best first) | aug-cc-pVTZ ranking (best first) |
  |---|---|---|
  | dipole | PBE (+0.030), B3LYP (+0.064), HF (+0.133) | B3LYP (-0.008), PBE (-0.050), HF (+0.084) |
  | angle | B3LYP (+0.02), PBE (-0.99), HF (+1.48) | PBE (-0.35), B3LYP (+0.58), HF (+1.80) |

  Without diffuse functions PBE looks best for the dipole, and with them B3LYP does. A single-basis comparison would have given the wrong ranking, so no method is declared "best" here. HF overestimates the dipole in both.

## Limitations

- This is **one molecule and one property set**, so it does not rank the methods in general. A proper benchmark would use many molecules.
- The experimental dipole is vibrationally averaged, while the computed values are for the equilibrium geometry, so differences of a few hundredths of a debye (for example B3LYP's -0.008 D) are not meaningful.
- The largest basis is aug-cc-pVTZ. No complete-basis-set extrapolation or quadruple-zeta check was done, so the basis-set limit of each property is not established.
- No dispersion corrections or solvent effects were included, and the DFT integration grid is the PySCF default.
- Total energies from different methods are not comparable, so only properties are compared across methods.

## Reproduce

PySCF runs on Linux, macOS and WSL, but **not natively on Windows**.

```bash
conda env create -f environment.yml
conda activate water-benchmark
python benchmark.py
```

The script writes tables (`basis_study.csv`, `method_study.csv`), a markdown summary and the two figures into `results/`. It takes about 70 seconds on a single CPU core. The committed results were generated with PySCF 2.15.1 and geomeTRIC 1.1.1; the version is also printed at the end of `results/summary.md`.

## Files

- `benchmark.py` runs both studies and writes everything in `results/`
- `environment.yml` lists the conda environment
- `results/` holds the generated tables and figures
- `LICENSE` gives the usage terms

## Usage & terms

This repository is shared for viewing and evaluation only. See `LICENSE`: no copying, modification, redistribution or reuse without written permission from the author. Everything is provided "as is", without warranty.

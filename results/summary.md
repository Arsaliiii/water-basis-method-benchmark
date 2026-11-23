## Part A: basis-set study (HF, fixed geometry)

| family | basis | basis_functions | energy_Eh | energy_drop_kcal_mol | dipole_D | dipole_error_D |
|---|---|---|---|---|---|---|
| Pople | sto-3g | 7 | -74.963063 |  | 1.726 | -0.129 |
| Pople | 6-31g | 13 | -75.983948 | 640.6 | 2.632 | +0.777 |
| Pople | 6-31g* | 19 | -76.010496 | 16.7 | 2.227 | +0.372 |
| Dunning | cc-pvdz | 24 | -76.026766 |  | 2.058 | +0.203 |
| Dunning | cc-pvtz | 58 | -76.057114 | 19.0 | 2.027 | +0.172 |
| aug-Dunning | aug-cc-pvdz | 41 | -76.041382 |  | 2.001 | +0.146 |
| aug-Dunning | aug-cc-pvtz | 92 | -76.060559 | 12.0 | 1.985 | +0.130 |

## Part B: method study (geometry optimised for each method and basis)

| basis | method | r_OH_A | r_OH_error_A | angle_deg | angle_error_deg | dipole_D | dipole_error_D |
|---|---|---|---|---|---|---|---|
| cc-pvtz | HF | 0.9406 | -0.0166 | 106.00 | +1.48 | 1.988 | +0.133 |
| cc-pvtz | PBE | 0.9697 | +0.0125 | 103.53 | -0.99 | 1.885 | +0.030 |
| cc-pvtz | B3LYP | 0.9613 | +0.0041 | 104.54 | +0.02 | 1.919 | +0.064 |
| aug-cc-pvtz | HF | 0.9410 | -0.0162 | 106.32 | +1.80 | 1.939 | +0.084 |
| aug-cc-pvtz | PBE | 0.9703 | +0.0131 | 104.17 | -0.35 | 1.805 | -0.050 |
| aug-cc-pvtz | B3LYP | 0.9618 | +0.0046 | 105.10 | +0.58 | 1.847 | -0.008 |

Experiment: r(OH) = 0.9572 A, angle = 104.52 deg, dipole = 1.855 D

Generated with PySCF 2.15.1

"""Water benchmark: how basis set and method choice change computed properties.

Part A  Basis-set study: Hartree-Fock at a fixed, near-experimental geometry,
        across Pople, Dunning and augmented Dunning basis sets.
Part B  Method study:    HF, PBE and B3LYP, each with its own optimised
                         geometry, at cc-pVTZ and at aug-cc-pVTZ.

Run:      python benchmark.py
Outputs:  results/  (CSV tables, summary.md, PNG figures)
Requires: PySCF, which runs on Linux, macOS or WSL (not native Windows).
"""

import tempfile
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # draw figures to files, no display needed

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyscf
from pyscf import dft, gto, scf
from pyscf.geomopt.geometric_solver import optimize

HARTREE_TO_KCAL = 627.509

# Starting geometry for water (O-H 0.958 A, H-O-H 104.4 deg, close to
# experiment). Part A uses it as is; Part B relaxes it for each method.
START_GEOMETRY = "O 0 0 0; H 0 0.757 0.587; H 0 -0.757 0.587"

# Part A basis sets, grouped by family. Energy drops are only reported inside
# a family, because only there is each basis a systematic improvement on the
# previous one.
BASIS_FAMILIES = {
    "Pople": ["sto-3g", "6-31g", "6-31g*"],
    "Dunning": ["cc-pvdz", "cc-pvtz"],
    "aug-Dunning": ["aug-cc-pvdz", "aug-cc-pvtz"],
}
# Pople sets are defined with 6 Cartesian d functions; PySCF defaults to 5
# spherical. Use Cartesian for Pople so results match the literature.
CARTESIAN_FAMILIES = {"Pople"}

# Part B: the same method comparison with and without diffuse functions.
OPT_BASES = ["cc-pvtz", "aug-cc-pvtz"]
METHODS = ["HF", "PBE", "B3LYP"]

# Experimental gas-phase water (NIST CCCBDB). The dipole is vibrationally
# averaged, so it is not exactly comparable to an equilibrium-geometry value.
EXPERIMENT = {"r_OH": 0.9572, "angle": 104.52, "dipole": 1.855}

RESULTS = Path("results")


def make_mf(mol, method):
    """Return an HF or DFT mean-field object for the given method name."""
    if method == "HF":
        return scf.RHF(mol)
    mf = dft.RKS(mol)
    mf.xc = method
    return mf


def dipole(mf):
    """Magnitude of the dipole moment in Debye."""
    return float(np.linalg.norm(mf.dip_moment(verbose=0)))


def geometry(mol):
    """Mean O-H bond length (Angstrom) and H-O-H angle (degrees)."""
    xyz = mol.atom_coords(unit="Angstrom")
    v1, v2 = xyz[1] - xyz[0], xyz[2] - xyz[0]
    n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
    angle = np.degrees(np.arccos(v1 @ v2 / (n1 * n2)))
    return 0.5 * (n1 + n2), float(angle)


def basis_study():
    """Part A: Hartree-Fock energy and dipole versus basis set."""
    rows = []
    for family, bases in BASIS_FAMILIES.items():
        for basis in bases:
            t0 = time.time()
            mol = gto.M(
                atom=START_GEOMETRY,
                basis=basis,
                cart=family in CARTESIAN_FAMILIES,
                verbose=0,
            )
            mf = scf.RHF(mol)
            energy = mf.kernel()
            rows.append(
                {
                    "family": family,
                    "basis": basis,
                    "basis_functions": mol.nao,
                    "energy_Eh": energy,
                    "dipole_D": dipole(mf),
                    "seconds": time.time() - t0,
                }
            )
    df = pd.DataFrame(rows)
    df["energy_drop_kcal_mol"] = (
        -df.groupby("family")["energy_Eh"].diff() * HARTREE_TO_KCAL
    )
    df["dipole_error_D"] = df["dipole_D"] - EXPERIMENT["dipole"]
    return df


# geomeTRIC prints a long banner for every optimisation. This log config only
# shows warnings so the console output stays readable.
QUIET_LOG = """[loggers]
keys=root

[handlers]
keys=stream_handler

[formatters]
keys=formatter

[logger_root]
level=WARNING
handlers=stream_handler

[handler_stream_handler]
class=geometric.nifty.RawStreamHandler
level=WARNING
formatter=formatter
args=(sys.stderr,)

[formatter_formatter]
format=%(message)s
"""


def method_study():
    """Part B: optimise the geometry with each method and basis, then compare."""
    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        log_ini = Path(tmp) / "quiet_log.ini"
        log_ini.write_text(QUIET_LOG)

        for basis in OPT_BASES:
            for method in METHODS:
                t0 = time.time()
                mol = gto.M(atom=START_GEOMETRY, basis=basis, verbose=0)
                mol_eq = optimize(make_mf(mol, method), logIni=str(log_ini))
                mf = make_mf(mol_eq, method)
                energy = mf.kernel()
                r_oh, angle = geometry(mol_eq)
                rows.append(
                    {
                        "basis": basis,
                        "method": method,
                        "energy_Eh": energy,
                        "r_OH_A": r_oh,
                        "angle_deg": angle,
                        "dipole_D": dipole(mf),
                        "seconds": time.time() - t0,
                    }
                )
    df = pd.DataFrame(rows)
    df["r_OH_error_A"] = df["r_OH_A"] - EXPERIMENT["r_OH"]
    df["angle_error_deg"] = df["angle_deg"] - EXPERIMENT["angle"]
    df["dipole_error_D"] = df["dipole_D"] - EXPERIMENT["dipole"]
    return df


def to_markdown(df, columns, formats):
    """Small markdown-table writer (avoids an extra dependency)."""
    lines = ["| " + " | ".join(columns) + " |", "|" + "---|" * len(columns)]
    for _, row in df.iterrows():
        cells = []
        for col in columns:
            value = row[col]
            if isinstance(value, str):
                cells.append(value)
            elif pd.isna(value):
                cells.append("")
            else:
                cells.append(format(value, formats.get(col, ".3f")))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def plot_basis(df):
    fig, ax = plt.subplots(figsize=(7, 4))
    colors = {"Pople": "#4c78a8", "Dunning": "#f58518", "aug-Dunning": "#54a24b"}
    for family, sub in df.groupby("family", sort=False):
        ax.plot(sub.index, sub["dipole_D"], "o-", color=colors[family], label=family)
    ax.axhline(EXPERIMENT["dipole"], color="k", ls="--", label="experiment")
    ax.set_xticks(df.index)
    ax.set_xticklabels(df["basis"], rotation=30, ha="right")
    ax.set_ylabel("Dipole moment (D)")
    ax.set_xlabel("Basis set")
    ax.set_title("Water dipole vs basis set (Hartree-Fock, fixed geometry)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(RESULTS / "dipole_vs_basis.png", dpi=150)
    plt.close(fig)


def plot_methods(df):
    panels = [
        ("r_OH_error_A", "O-H bond error (A)"),
        ("angle_error_deg", "H-O-H angle error (deg)"),
        ("dipole_error_D", "Dipole error (D)"),
    ]
    colors = {"cc-pvtz": "#999999", "aug-cc-pvtz": "#4c78a8"}
    width = 0.38
    x = np.arange(len(METHODS))
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.8))
    for ax, (col, label) in zip(axes, panels):
        for i, basis in enumerate(OPT_BASES):
            sub = df[df["basis"] == basis].set_index("method").loc[METHODS]
            ax.bar(x + (i - 0.5) * width, sub[col], width,
                   color=colors[basis], label=basis)
        ax.axhline(0, color="k", lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(METHODS)
        ax.set_title(label, fontsize=10)
    axes[0].legend(fontsize=8)
    fig.suptitle("Error vs experiment, geometry optimised with each method and basis")
    fig.tight_layout()
    fig.savefig(RESULTS / "method_errors.png", dpi=150)
    plt.close(fig)


def main():
    RESULTS.mkdir(exist_ok=True)

    print("Part A: basis-set study (HF, fixed geometry)")
    basis_df = basis_study()
    basis_df.to_csv(RESULTS / "basis_study.csv", index=False)

    print("Part B: method study (geometry optimised for each method and basis)")
    method_df = method_study()
    method_df.to_csv(RESULTS / "method_study.csv", index=False)

    plot_basis(basis_df)
    plot_methods(method_df)

    basis_md = to_markdown(
        basis_df,
        ["family", "basis", "basis_functions", "energy_Eh",
         "energy_drop_kcal_mol", "dipole_D", "dipole_error_D"],
        {"basis_functions": ".0f", "energy_Eh": ".6f",
         "energy_drop_kcal_mol": ".1f", "dipole_error_D": "+.3f"},
    )
    method_md = to_markdown(
        method_df,
        ["basis", "method", "r_OH_A", "r_OH_error_A", "angle_deg",
         "angle_error_deg", "dipole_D", "dipole_error_D"],
        {"r_OH_A": ".4f", "r_OH_error_A": "+.4f", "angle_deg": ".2f",
         "angle_error_deg": "+.2f", "dipole_error_D": "+.3f"},
    )
    summary = (
        "## Part A: basis-set study (HF, fixed geometry)\n\n"
        f"{basis_md}\n\n"
        "## Part B: method study (geometry optimised for each method and basis)\n\n"
        f"{method_md}\n\n"
        f"Experiment: r(OH) = {EXPERIMENT['r_OH']} A, "
        f"angle = {EXPERIMENT['angle']} deg, dipole = {EXPERIMENT['dipole']} D\n\n"
        f"Generated with PySCF {pyscf.__version__}\n"
    )
    (RESULTS / "summary.md").write_text(summary)
    print("\n" + summary)
    print(f"Saved tables and figures in {RESULTS}/")


if __name__ == "__main__":
    main()

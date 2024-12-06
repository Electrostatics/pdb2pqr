"""Common routines and variables for testing"""

import hashlib
import itertools
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from pdb2pqr.main import build_main_parser, main_driver
from pdb2pqr.structures import Atom

_LOGGER = logging.getLogger(__name__)
PARSER = build_main_parser()


# This is a list of legit FF names
FF_LIST = ["AMBER", "CHARMM", "PARSE", "TYL06", "PEOEPB", "SWANSON"]
# Where the test data lives
DATA_DIR = Path("tests/data")
# Tolerable error for positions
POS_CUT = 1e-2
# Tolerable error for charge
Q_CUT = 1e-2
# Tolerable error for radius
R_CUT = 1e-2


def generate_ff_combinations():
    """Generate combinations of force field options.

    Returns:
        list of combinations
    """
    ff_options = [None] + [f"--ff={ff}" for ff in FF_LIST]
    ffout_options = [None] + [f"--ffout={ff}" for ff in FF_LIST]
    ff_combos = []
    for ff_opt in ff_options:
        for ffout_opt in ffout_options:
            combo = []
            if ff_opt is not None:
                combo += [ff_opt]
            if ffout_opt is not None:
                combo += [ffout_opt]
            ff_combos += [combo]
    return ff_combos


def generate_combinations(option_list):
    """Generate all combinations of provided options.

    Returns:
        List of option combinations.
    """
    combos = []
    for num in range(len(option_list) + 1):
        combos += list(itertools.combinations(option_list, num))
    return combos


def pqr_to_dict(pqr_file):
    """Convert PQR to dictionary.

    Args:
        pqr_file:  file object for PQR.

    Returns:
        DataFrame with PQR information.
    """
    pqr = []
    for line in pqr_file:
        atom = Atom.from_pqr_line(line)
        if atom is not None:
            row_dict = {"atom_num": atom.serial}
            # Many hydrogens are created in arbitrary order when attached to
            # the same heavy atom. Therefore, the last number in their name is
            # not meaningful
            atom_name = atom.name
            if atom_name[0] == "H":
                try:
                    int(atom_name[-1])
                    atom_name = atom_name[:-1] + "_"
                except ValueError:
                    pass
            row_dict["atom_name"] = atom_name
            row_dict["res_name"] = atom.res_name
            row_dict["res_num"] = atom.res_seq
            row_dict["x"] = atom.x
            row_dict["y"] = atom.y
            row_dict["z"] = atom.z
            row_dict["q"] = atom.charge
            row_dict["r"] = atom.radius
            pqr.append(row_dict)
    return pd.DataFrame(pqr)


def pqr_distance(df1, df2):
    """Calculate distances between positions, charges, and radii from two PQR dataframes.

    Args:
        df1:  PQR dataframe
        df2:  PQR dataframe

    Returns:
        Dataframe of distances
    """
    if "chain" in df1.columns:
        d_frame = df1.merge(
            df2,
            on=["atom_name", "res_name", "res_num", "chain"],
            how="inner",
            suffixes=("A", "B"),
        )
    else:
        d_frame = df1.merge(
            df2,
            on=["atom_name", "res_name", "res_num"],
            how="inner",
            suffixes=("A", "B"),
        )

    # Calculate differences and drop original columns
    for c_val in ("x", "y", "z", "q", "r"):
        d_val = f"d{c_val}"
        d2_val = f"d{c_val}2"
        c_a = f"{c_val}A"
        c_b = f"{c_val}B"
        d_frame[d_val] = d_frame[c_a] - d_frame[c_b]
        d_frame[d2_val] = d_frame[d_val] * d_frame[d_val]
        d_frame = d_frame.drop([d_val, c_a, c_b], axis="columns")

    # Calculate position norm-squared and drop used columns
    d_frame["dp2"] = d_frame["dx2"] + d_frame["dy2"] + d_frame["dz2"]
    d_frame = d_frame.drop(["dx2", "dy2", "dz2"], axis="columns")

    # Calculate norms of all measures and drop used columns
    for c_val in ("p", "q", "r"):
        n_val = f"d{c_val}"
        n2_val = f"d{c_val}2"
        d_frame[n_val] = np.sqrt(d_frame[n2_val])
        d_frame = d_frame.drop(n2_val, axis="columns")

    return d_frame


def compare_pqr(pqr1_path, pqr2_path, compare_resnames=False):
    """Compare two PQR files.

    Assume that atom numbering/ordering does not matter.

    Args:
        pqr1_path:  Path to first PQR
        par2_path:  Path to second PQR
    """
    with open(pqr1_path, encoding="utf-8") as pqr1_file:
        df1 = pqr_to_dict(pqr1_file)
        _LOGGER.debug(f"PQR 1 has shape {df1.shape}")

    with open(pqr2_path, encoding="utf-8") as pqr2_file:
        df2 = pqr_to_dict(pqr2_file)
        _LOGGER.debug(f"PQR 2 has shape {df2.shape}")

    if compare_resnames:
        diff_resname = df1.res_name != df2.res_name
        if any(diff_resname):
            _LOGGER.error(f"Got residue names:\n{df1[diff_resname]}")
            _LOGGER.error(f"Expected residue names:\n{df2[diff_resname]}")
            raise ValueError(
                f"PQR files have different residue names\n{df1.res_name[diff_resname]}\nreference\n{df2.res_name[diff_resname]}"
            )

    d_frame = pqr_distance(df1, df2)
    _LOGGER.debug(f"Merged d_frame has shape {d_frame.shape}")

    grouped = d_frame.groupby(["res_name", "res_name", "atom_name"])
    _LOGGER.debug(f"Have {len(grouped):d} unique atoms")
    df_min = grouped.min()

    for col, what, cut in [
        ("dp", "position", POS_CUT),
        ("dq", "charge", Q_CUT),
        ("dr", "radius", R_CUT),
    ]:
        for cut_ in [0.0, cut]:
            df_c = df_min[df_min[col] > cut_].sort_values(col, ascending=False)
            ndiff = df_c.shape[0]
            result = f"{ndiff:d} atoms have {what} differences > {cut_:g}"
            if ndiff > 0:
                _LOGGER.warning(result)
                df_c = df_min[df_min[col] > cut_].sort_values(
                    col, ascending=False
                )
                summary = [
                    f"{key}: {val:.3E}"
                    for (key, val) in df_c[col].describe().to_dict().items()
                ]
                _LOGGER.debug(summary)
                if cut_ > 0:
                    raise ValueError(result)
            else:
                _LOGGER.info(result)


def run_pdb2pqr_for_tests(
    args,
    input_pdb,
    tmp_path,
    output_pqr=None,
    expected_pqr=None,
    compare_resnames=False,
):
    """Basic code for invoking PDB2PQR."""
    from pdb2pqr import io

    if output_pqr is None:
        hash_str = f"{args}{input_pdb}"
        hash_ = hashlib.sha1(hash_str.encode("UTF-8")).hexdigest()
        output_pqr = f"{hash_}.pqr"
    output_pqr = tmp_path / output_pqr
    _LOGGER.debug(f"Writing output to {output_pqr}")
    arg_str = f"{args} {input_pdb} {output_pqr}"
    args = PARSER.parse_args(arg_str.split())
    io.setup_logger(args.output_pqr, args.log_level)
    _ = main_driver(args)
    if expected_pqr is not None:
        compare_pqr(
            output_pqr, expected_pqr, compare_resnames=compare_resnames
        )


def run_propka_for_tests(input_pdb, compare_file, pH):
    from propka.lib import build_parser

    from pdb2pqr import io
    from pdb2pqr.main import run_propka, setup_molecule

    propka_parser = build_parser()
    propka_parser.add_argument(
        "--keep-chain", action="store_true", default=False
    )
    propka_args = propka_parser.parse_args(
        ["--log-level", "WARNING", "--pH", str(pH), "TEST"]
    )

    definition = io.get_definitions()
    pdblist, _ = io.get_molecule(input_pdb)
    biomolecule, definition, _ = setup_molecule(pdblist, definition, None)
    biomolecule.set_termini(neutraln=False, neutralc=False)
    biomolecule.update_bonds()

    biomolecule.remove_hydrogens()
    pka_df, _ = run_propka(propka_args, biomolecule)
    results = {
        f"{row['res_name']} {row['res_num']} {row['chain_id']}": row["pKa"]
        for row in pka_df
        if row["group_label"].startswith(row["res_name"])
    }

    # # Used for writing the test CSV files
    # with open(compare_file, "w") as f:
    #     f.write("group,pka\n")
    #     for key, val in results.items():
    #         f.write(f"{key},{val}\n")

    compare = {}
    with open(compare_file) as f:
        for line in f.readlines()[1:]:
            line = line.strip()
            if len(line) == 0:
                continue
            key, val = line.split(",")
            compare[key] = float(val)

    diff_new = [key for key in results if key not in compare]
    if len(diff_new):
        raise RuntimeError(
            f"New pKa groups {' '.join(diff_new)} were generated by test which didn't exist in the test reference files."
        )
    diff_old = [key for key in compare if key not in results]
    if len(diff_old):
        raise RuntimeError(
            f"Missing pKa groups {' '.join(diff_old)} which exist in the test reference files."
        )

    for key in results:
        if abs(results[key] - compare[key]) > 1e-3:
            raise RuntimeError(
                f"PROPKA test error. pKa values changed more than 1e-3 for group {key}: new {results[key]} reference {compare[key]}"
            )

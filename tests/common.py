"""Common routines and variables for testing"""
import logging
import itertools
import hashlib
from pathlib import Path
import numpy
import pandas
from pdb2pqr.main import build_parser, main


_LOGGER = logging.getLogger(__name__)
PARSER = build_parser()


# This is a list of legit FF names
FF_LIST = ["AMBER", "CHARMM", "PARSE", "TYL06", "PEOPB", "SWANSON"]
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
        list of combinations"""
    ff_options = [None] + ["--ff=%s" % ff for ff in FF_LIST]
    ffout_options = [None] + ["--ffout=%s" % ff for ff in FF_LIST]
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
    for num in range(len(option_list)+1):
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
        row_dict = {}
        line = " ".join([line[:6], line[6:]])
        words = line.strip().split()
        label = words.pop(0)
        if label in ["REMARK", "TER", "END", "HEADER", "TITLE", "COMPND",
                     "SOURCE", "KEYWDS", "EXPDTA", "AUTHOR", "REVDAT", "JRNL"]:
            pass
        elif label in ["ATOM", "HETATM"]:
            row_dict["atom_num"] = int(words.pop(0))
            # Many hydrogens are created in arbitrary order when attached to
            # the same heavy atom. Therefore, the last number in their name is
            # not meaningful
            atom_name = words.pop(0).strip()
            if atom_name[0] == "H":
                try:
                    int(atom_name[-1])
                    atom_name = atom_name[:-1] + "_"
                except ValueError:
                    pass
            row_dict["atom_name"] = atom_name
            row_dict["res_name"] = words.pop(0).strip()
            num_or_chain = words.pop(0).strip()
            try:
                res_num = int(num_or_chain)
                row_dict["res_num"] = res_num
            except ValueError:
                row_dict["chain"] = num_or_chain
                row_dict["res_num"] = int(words.pop(0))
            row_dict["x"] = float(words.pop(0))
            row_dict["y"] = float(words.pop(0))
            row_dict["z"] = float(words.pop(0))
            row_dict["q"] = float(words.pop(0))
            row_dict["r"] = float(words.pop(0))
            pqr.append(row_dict)
        else:
            raise NotImplementedError(label, words)
    return pandas.DataFrame(pqr)


def pqr_distance(df1, df2):
    """Calculate distances between positions, charges, and radii from two PQR
    dataframes.

    Args:
        df1:  PQR dataframe
        df2:  PQR dataframe
    
    Returns:
        Dataframe of distances
    """
    if "chain" in df1.columns:
        df = df1.merge(df2, on=["atom_name", "res_name", "res_num", "chain"],
                       how="inner", suffixes=("A", "B"))
    else:
        df = df1.merge(df2, on=["atom_name", "res_name", "res_num"], how="inner",
                       suffixes=("A", "B"))

    # Calculate differences and drop original columns
    for c in ("x", "y", "z", "q", "r"):
        d = "d%s" % c
        d2 = "d%s2" % c
        cA = "%sA" % c
        cB = "%sB" % c
        df[d] = df[cA] - df[cB]
        df[d2] = df[d] * df[d]
        df = df.drop([d, cA, cB], axis="columns")

    # Calculate position norm-squared and drop used columns
    df["dp2"] = df["dx2"] + df["dy2"] + df["dz2"]
    df = df.drop(["dx2", "dy2", "dz2"], axis="columns")

    # Calculate norms of all measures and drop used columns
    for c in ("p", "q", "r"):
        n = "d%s" % c
        n2 = "d%s2" % c
        df[n] = numpy.sqrt(df[n2])
        df = df.drop(n2, axis="columns")
    
    return df


def compare_pqr(pqr1_path, pqr2_path):
    """Compare two PQR files.

    Assume that atom numbering/ordering does not matter.

    Args:
        pqr1_path:  Path to first PQR
        par2_path:  Path to second PQR
    """
    with open(pqr1_path, "rt", encoding="utf-8") as pqr1_file:
        df1 = pqr_to_dict(pqr1_file)
        _LOGGER.debug("PQR 1 has shape %s", df1.shape)
    
    with open(pqr2_path, "rt", encoding="utf-8") as pqr2_file:
        df2 = pqr_to_dict(pqr2_file)
        _LOGGER.debug("PQR 2 has shape %s", df2.shape)
    
    df = pqr_distance(df1, df2)
    _LOGGER.debug("Merged df has shape %s", df.shape)

    grouped = df.groupby(["res_name", "res_name", "atom_name"])
    _LOGGER.debug("Have %d unique atoms", len(grouped))
    df_min = grouped.min()

    for col, what, cut in [("dp", "position", POS_CUT),
                           ("dq", "charge", Q_CUT),
                           ("dr", "radius", R_CUT)]:
        for cut_ in [0.0, cut]:
            df_c = df_min[df_min[col] > cut_].sort_values(col, ascending=False)
            ndiff = df_c.shape[0]
            result = "%d atoms have %s differences > %g" % (ndiff, what, cut_)
            if ndiff > 0:
                _LOGGER.warning(result)
                df_c = df_min[df_min[col] > cut_].sort_values(col, ascending=False)
                summary = ["%s: %.3E" % (key, val) for (key, val) in df_c[col].describe().to_dict().items()]
                _LOGGER.debug(summary)
                if cut_ > 0:
                    raise ValueError(result)
            else:
                _LOGGER.info(result)


def run_pdb2pqr(args, input_pdb, tmp_path, output_pqr=None, expected_pqr=None):
    """Basic code for invoking PDB2PQR."""
    arg_str = args + " {inp} {out}"
    if output_pqr is None:
        hash_str = str(args) + str(input_pdb)
        hash_ = hashlib.sha1(hash_str.encode("UTF-8")).hexdigest()
        output_pqr = hash_ + ".pqr"
    output_pqr = tmp_path / output_pqr
    _LOGGER.debug("Writing output to %s", output_pqr)
    arg_str = arg_str.format(inp=input_pdb, out=output_pqr)
    args = PARSER.parse_args(arg_str.split())
    main(args)
    if expected_pqr is not None:
        compare_pqr(output_pqr, expected_pqr)

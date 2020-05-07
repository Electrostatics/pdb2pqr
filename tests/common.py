"""Common routines and variables for testing"""
from pathlib import Path
import numpy
import pandas
import logging


_LOGGER = logging.getLogger(__name__)


DATA_DIR = Path("tests/data")
POS_CUT = 1e-2
Q_CUT = 1e-2
R_CUT = 1e-2

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
                _LOGGER.info(summary)
                if cut_ > 0:
                    raise ValueError(result)
            else:
                _LOGGER.info(result)
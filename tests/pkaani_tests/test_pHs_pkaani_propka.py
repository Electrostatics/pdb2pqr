import sys
import subprocess
from pathlib import Path
from collections import Counter


def count_residues(file_path):
    """
    Count residue types by reading lines directly.
    Only considers ATOM lines with alpha carbon (CA) atoms.
    """
    counts = Counter()
    with open(file_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                residue = line[17:20].strip()
                counts[residue] += 1
    return counts


def run_pdb2pqr(pdb_file, method, pH):
    base_name = Path(pdb_file).stem
    output_file = f"{base_name}_{method}.pqr"
    subprocess.run([
        "pdb2pqr30", "--ffout", "AMBER",
        "--titration-state-method", method,
        "--with-ph", pH,
        str(pdb_file), output_file
    ], check=True)
    return output_file


def print_results(title, original_counts, pqr_counts):
    def safe_get(d, key):
        return d.get(key, 0)

    print(f"{title}")
    print("---------------------------------------")
    print("Aspartates in pdb file: (most should not be protonated after pdb2pqr at pH 5.0)")
    print(safe_get(original_counts, "ASP"))
    print("Protonated aspartates in pqr file:")
    print(safe_get(pqr_counts, "ASH"))

    print("Glutamates in pdb file: (most should not be protonated after pdb2pqr at pH 5.0)")
    print(safe_get(original_counts, "GLU"))
    print("Protonated glutamates in pqr file:")
    print(safe_get(pqr_counts, "GLH"))

    print("Histidines in pdb file: (all should be protonated after pdb2pqr at pH 5.0)")
    print(safe_get(original_counts, "HIS"))
    print("Protonated histidines in pqr file:")
    print(safe_get(pqr_counts, "HIP"))

    print("Protonated Lysines in pdb file:")
    print(safe_get(original_counts, "LYS"))
    print("Lysines in pqr file:")
    print(safe_get(pqr_counts, "LYS"))

    print("Protonated Arginines in pdb file:")
    print(safe_get(original_counts, "ARG"))
    print("Arginines in pqr file:")
    print(safe_get(pqr_counts, "ARG"))

    print("Protonated Tyrosines in pdb file:")
    print(safe_get(original_counts, "TYR"))
    print("Tyrosines in pqr file:")
    print(safe_get(pqr_counts, "TYR"))
    print("---------------------------------------\n")


def main():
    if len(sys.argv) != 3:
        print("Usage: python test_pHs_pkaani_propka.py input.pdb pH")
        sys.exit(1)

    pdb_file = Path(sys.argv[1])
    pH = sys.argv[2]

    if not pdb_file.is_file():
        print(f"File not found: {pdb_file}")
        sys.exit(1)

    # Count residues in original PDB
    original_counts = count_residues(pdb_file)

    for method in ["propka", "pkaani"]:
        pqr_file = run_pdb2pqr(pdb_file, method, pH)
        pqr_counts = count_residues(pqr_file)
        print_results(f"{method.upper()} RESULTS {pdb_file.name}", original_counts, pqr_counts)


if __name__ == "__main__":
    main()
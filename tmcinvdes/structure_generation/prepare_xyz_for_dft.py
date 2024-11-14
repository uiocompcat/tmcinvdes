"""Prepare XYZs of assembled TMCs into individual subfolders and files for DFT
labeling.

python -m tmcinvdes.structure_generation.prepare_xyz_for_dft \
                -i datasets/12_cond-TMC/cond_mono-isolated_ligands-sampled_optimized-TMC.xyz \
                -o tmp_results/orca_xyz-cond_mono-isolated_ligands-sampled_optimized-TMC
"""

import argparse
import os
from pathlib import Path


def parse_args(arg_list: list = None) -> argparse.Namespace:
    """Parse arguments from command line.

    Args:
        arg_list (list, optional): Automatically obtained from the command line if provided.
        If no arguments are given but default arguments are defined, the latter are used.

    Returns:
        argparse.Namespace: Dictionary-like class that contains the driver arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_file",
        "-i",
        type=Path,
        required=True,
        help="Input file with multiple XYZ blocks.",
    )
    parser.add_argument(
        "--output_dir",
        "-o",
        type=Path,
        required=True,
        help="Directory of output files.",
    )
    return parser.parse_args(arg_list)


if __name__ == "__main__":
    args = parse_args()
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=False)
    input_file = os.path.abspath(args.input_file)

    xyzs = []
    with open(input_file, "r") as f:
        xyzs = f.read().split("\n\n")
    for i, xyz in enumerate(xyzs):
        xyz_lines = [line.rstrip() for line in xyz.rstrip().split("\n") if line != ""]
        filename = xyz_lines[1] + ".xyz"
        sample_dir = os.path.join(output_dir, f"sample-{i}")
        os.makedirs(sample_dir, exist_ok=False)
        filepath = os.path.join(sample_dir, filename)
        print(filepath)
        with open(filepath, "w") as f:
            f.write("\n".join(xyz_lines))

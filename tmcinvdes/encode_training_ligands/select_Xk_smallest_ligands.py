"""Select a subset of smaller ligands and give them IDs as part of a new set of
ligands.

Usage:

python select_Xk_smallest_ligands.py -d denticity -i input_dir -o output_dir -x test

python select_Xk_smallest_ligands.py -d bidentate \
                                     -i ../../datasets/04_uncond_novel \
                                     -o ../../datasets/05_uncond_minXk \
                                     -x test
"""

import argparse
import os
from pathlib import Path

import pandas as pd

from tmcinvdes.ligand_generation.utils import compare_dataframes


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
        "--denticity",
        "-d",
        choices=[
            "monodentate",
            "bidentate",
        ],
        default="bidentate",
        help="Select one of the two denticity modes supported.",
    )
    parser.add_argument(
        "--input_dir",
        "-i",
        type=Path,
        required=True,
        help="Input directory with the initial screen of the generated ligands.",
    )
    parser.add_argument(
        "--output_dir",
        "-o",
        type=Path,
        required=True,
        help="Output directory for the selected subset of Xk smallest ligands",
    )
    parser.add_argument(
        "--xtent",
        "-x",
        choices=["full", "test"],
        default="test",
        help="""The extent to which the process should run.
         - In the "full" case, the process runs in its entirety and writes the output results to
           file.
         - In the "test" case, the process runs in its entirety and compares the results to the
           output already existing in the designated output file.""",
    )
    return parser.parse_args(arg_list)


def select_ligands_subset(
    df_input: pd.DataFrame, xk: int, dataset_name: str
) -> pd.DataFrame:
    """Select the first X thousand ligands.

    Args:
        df_input (pd.DataFrame): dataframe of ligands already ordered by (heavy) atom counts,
        depending on denticity.
        xk (int): the X thousand target number of ligands to select.

    Returns:
        pd.DataFrame: the selected subset of ligands.
    """
    df_output = df_input[:xk]
    ligand_ids = []
    for i in range(df_output.shape[0]):
        ligand_id = f"{dataset_name}-{i+1}"
        ligand_ids.append(ligand_id)
    df_output.insert(0, "Ligand ID", ligand_ids)
    return df_output


if __name__ == "__main__":
    args = parse_args()
    denticity = args.denticity
    xtent = args.xtent
    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output_dir)

    # Modify parameters based on denticity:
    if denticity == "monodentate":
        input_file = os.path.join(input_dir, "uncond_mono-raw50k-novel.csv")
        output_file = os.path.join(output_dir, "uncond_mono-min15k.csv")
        atom_counting_mode = "Heavy atom count"
        xk = 15000
        dataset_name = "uncond_mono-min15k"
    elif denticity == "bidentate":
        input_file = os.path.join(input_dir, "uncond_bi-raw50k-novel.csv")
        output_file = os.path.join(output_dir, "uncond_bi-min10k.csv")
        atom_counting_mode = "Atom count"
        xk = 10000
        dataset_name = "uncond_bi-min10k"

    df_input = pd.read_csv(input_file)

    # Modify process based on xtent:
    if xtent == "test":
        df_expect = pd.read_csv(output_file)
        df_output = select_ligands_subset(df_input, xk, dataset_name)
        row_accuracy = compare_dataframes(df_output, df_expect, atom_counting_mode)
        if row_accuracy == 1.0:
            print(
                f"Test run perfectly reproduced the pre-existing file: {output_file}."
            )
        else:
            print(
                f"Test run failed to perfectly reproduce the pre-existing file: {output_file}."
            )
            print(
                f"Reproduction accuracy by matching rows over total rows: {row_accuracy}."
            )
    elif xtent == "full":
        df_output = select_ligands_subset(df_input, xk, dataset_name)
        df_output.to_csv(output_file, index=False, header=True)

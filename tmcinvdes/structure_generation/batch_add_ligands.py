"""Add a batch of ligands to the local molSimplify ligand library.

Usage:(from root directory of repository)

python -m tmcinvdes.structure_generation.batch_add_ligands -d bidentate \
                                                -i datasets/05_uncond_minXk/uncond_bi-min10k.csv \
                                                -s 0
"""

import argparse
from pathlib import Path

import pandas as pd
from rdkit import Chem

from tmcinvdes.ligand_generation.utils import (
    process_substitute_attachment_points,
    process_substitute_attachment_points_bidentate,
)
from tmcinvdes.structure_generation.molsimplify_tools import (
    add_ligand_from_smiles,
    get_existing_ligand_names,
)


def parse_args(arg_list: list = None) -> argparse.Namespace:
    """Parse arguments from command line.

    Args:
        arg_list (list, optional): Automatically obtained from the command line if provided.
        If no arguments are given but default arguments are defined, the latter are used.

    Returns:
        argparse.Namespace: Dictionary-like class that contains the driver arguments.
    """
    parser = argparse.ArgumentParser(
        description="Strip substitute atoms from generated ligands"
    )
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
        "--input_file",
        "-i",
        type=Path,
        default="datasets/uncond_bi_min10k/uncond_bi_min10k.csv",
        required=True,
        help="Input file with ligand data, with a minimum column set required.",
    )
    parser.add_argument(
        "--optimized",
        "-p",
        type=bool,
        default=False,
        help="Indicate whether ligands in input file were optimized.",
    )
    parser.add_argument(
        "--start",
        "-s",
        type=int,
        default=0,
        required=True,
        help="Offset starting point when iterating over ligands in input file.",
    )
    return parser.parse_args(arg_list)


def batch_add_ligands(
    df_ligands: pd.DataFrame, start: int, denticity: str = "monodentate"
):
    """Simplify the ligand addtion process by checking for pre-existing ligands
    only once per batch.

    Args:
        df_ligands (pd.DataFrame): the input dataframe corresponding to stage 6 with ligands to
        add.
        start (int): offset index if the input dataframe was already partly processed.
        denticity (str): the denticity of all the ligands in the input dataframe.
    """
    extant_ligands = set(get_existing_ligand_names())
    pre_checked = True
    for i, ligand_id in enumerate(df_ligands["Ligand ID"].values):
        if i < start:
            continue
        if ligand_id in extant_ligands:
            continue
        print(f"{i}: {ligand_id}")
        df_temp = df_ligands[df_ligands["Ligand ID"] == ligand_id]
        if "Decoded SMILES" in df_temp.columns:
            decoded_smiles = df_temp["Decoded SMILES"].values[0]
            connection_ids = eval(df_temp["Connection IDs"].values[0])
        else:
            encoded_smiles = df_temp["Encoded SMILES"].values[0]
            encoded_mol = Chem.MolFromSmiles(encoded_smiles)
            if denticity == "monodentate":
                decoded_mol, connection_id = process_substitute_attachment_points(
                    encoded_mol
                )
                connection_ids = [connection_id]
            elif denticity == "bidentate":
                (
                    decoded_mol,
                    connection_ids,
                ) = process_substitute_attachment_points_bidentate(encoded_mol)
            if decoded_mol is None:
                print(
                    f"Failed to process substitute attachment points of encoded {ligand_id}."
                )
                print(encoded_smiles)
                print("")
                continue
            else:
                decoded_smiles = Chem.MolToSmiles(decoded_mol)

        add_ligand_from_smiles(
            decoded_smiles, ligand_id, connection_ids, pre_checked=pre_checked
        )


def batch_add_optimized_ligands(df_ligands: pd.DataFrame, start: int, denticity: str):
    """Simplify the ligand addtion process by checking for pre-existing ligands
    only once per batch.

    This variant function takes care of differing column names used in stage 11 of the workflow.

    Args:
        df_ligands (pd.DataFrame): the input dataframe corresponding to stage 11 with ligands to
        add.
        start (int): offset index if the input dataframe was already partly processed.
        denticity (str): the denticity of all the optimized ligands.
    """
    extant_ligands = set(get_existing_ligand_names())
    pre_checked = True
    for i, ligand_id in enumerate(df_ligands["Optimized ligand ID"].values):
        if i < start:
            continue
        if ligand_id in extant_ligands:
            continue
        print(f"{i}: {ligand_id}")
        valid = False
        df_temp = df_ligands[df_ligands["Optimized ligand ID"] == ligand_id]
        encoded_smiles = df_temp["Optimized encoded SMILES"].values[0]
        encoded_mol = Chem.MolFromSmiles(encoded_smiles)
        if denticity == "monodentate":
            decoded_mol, connection_id = process_substitute_attachment_points(
                encoded_mol
            )
            connection_ids = [connection_id]
        elif denticity == "bidentate":
            (
                decoded_mol,
                connection_ids,
            ) = process_substitute_attachment_points_bidentate(encoded_mol)
        if decoded_mol is None:
            valid = False
            print(
                f"Failed to process substitute attachment points of encoded {ligand_id}."
            )
            print(encoded_smiles)
            print("")
            continue
        else:
            decoded_smiles = Chem.MolToSmiles(decoded_mol)
            valid = True
        if valid:
            add_ligand_from_smiles(
                decoded_smiles, ligand_id, connection_ids, pre_checked=pre_checked
            )


if __name__ == "__main__":
    args = parse_args()
    ligands_file = str(args.input_file)
    df_ligands = pd.read_csv(ligands_file)
    start = int(args.start)
    denticity = args.denticity

    if args.optimized:
        batch_add_optimized_ligands(df_ligands, start, denticity)
    else:
        batch_add_ligands(df_ligands, start, denticity)

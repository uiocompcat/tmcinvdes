"""Script for obtaining SMILES representations of tmQMg-L ligands used for
training the JT-VAE."""

import argparse
import multiprocessing as mp
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd
import swifter  # ruff: noqa: F401
import xyz2mol_functionality as x2m_tm
from rdkit import Chem
from utils import (
    attach_dummy_atom_to_coordinating_atoms,
    get_bidentate,
    get_connection_ids,
    get_monodentate,
    get_stable_occ,
    load_ligand_xyz,
    prune_num_atoms,
)

# Current file location to ensure paths do not break.
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


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
        default="monodentate",
        help="Select one of the two denticity modes supported.",
    )
    parser.add_argument(
        "--input_dir",
        "-i",
        type=Path,
        required=True,
        help="Input directory with CSV files from tmQMg-L",
    )
    parser.add_argument(
        "--output_dir",
        "-o",
        type=Path,
        required=True,
        default=".",
        help="""Output directory for training CSV file""",
    )
    parser.add_argument("--debug", action="store_true")

    return parser.parse_args(arg_list)


def ligand_xyz_to_mol(row: tuple, df_stable):
    """Function used to obtain ligand mol objects. The input is a row from a
    dataframe of the tmQMg-L dataset. It was used as a pandas apply function.

    Args:
        row (tuple): a row of a Pandas dataframe constructed from tmQMg-L CSV files.

    Returns:
        mol : RDKit mol object
        connect_ids : coordination atom idx
        stable_occ : Label of most stable ligand xyz file
    """
    row = row[1]
    connect_ids = get_connection_ids(row, df_stable)

    stable_occ = get_stable_occ(row["name"], df_stable)

    xyz = ligand_xyzs[stable_occ]

    # There is a bug with the charge column that changes its type
    if isinstance(row["charge"], np.integer):
        charge = row["charge"].item()
    else:
        charge = row["charge"]

    # Check if we can get mol object from tmQMg-L SMILES
    try:
        m = Chem.MolFromSmiles(row["smiles"])
        if not m:
            # If not, we assume the ligand is an edge case and discard
            return None
    except Exception:
        return None

    # Check for radical. If it contains radicals, we probably need xyz2mol.
    flag = False
    for atom in m.GetAtoms():
        if atom.GetNumRadicalElectrons() == 1:
            flag = True

    if flag:
        mol = x2m_tm.get_mol(xyz, charge=charge)

    # Check for carbenes.
    elif m.GetSubstructMatch(
        Chem.MolFromSmarts("[#6&v2H0,#6&v3H0]")
    ) or m.GetSubstructMatch(Chem.MolFromSmarts("[#14&v2H0]")):
        print("babel")
        mol = x2m_tm.get_mol_pure_babel(xyz, charge=charge)
    # Ligands with these elements get better SMILES with OpenBabel
    elif any(x in row["smiles"] for x in ["As", "Se", "p", "P"]):
        mol = x2m_tm.get_mol_pure_babel(xyz, charge=charge)
    else:
        mol = x2m_tm.get_mol(xyz, charge=charge)
    if not mol:
        return None

    # If there are still radicals present discard ligand
    for atom in mol.GetAtoms():
        if atom.GetNumRadicalElectrons() == 1:
            return None

    # Check the valency of coordinating atoms
    for connect_id in connect_ids:
        connect_atom = mol.GetAtomWithIdx(connect_id[0])
        id = connect_atom.GetIdx()

        valence = mol.GetAtomWithIdx(id).GetExplicitValence()

        # P with 5 bonds is assumed to be an edge case and discarded
        if connect_atom.GetSymbol() == "P":
            if valence > 4:
                return None

        # Prevent pentavelent coordinating carbons.
        if connect_atom.GetSymbol() == "C":
            if valence > 3:
                return None
            # Anionic carbons are discarded
            if connect_atom.GetFormalCharge() == -1:
                return None

    return (mol, connect_ids, stable_occ)


if __name__ == "__main__":
    args = parse_args()
    denticity = args.denticity

    args.output_dir.mkdir(exist_ok=True, parents=True)

    input_dir = os.path.abspath(args.input_dir)

    # Load tmQMg-L data
    df_desc = pd.read_csv(os.path.join(input_dir, "ligands_descriptors.csv"), sep=";")
    df_fingerprints = pd.read_csv(
        os.path.join(input_dir, "ligands_fingerprints.csv"), sep=";"
    )
    df_misc = pd.read_csv(os.path.join(input_dir, "ligands_misc_info.csv"), sep=";")
    df_stable = pd.read_csv(os.path.join(input_dir, "stable.csv"), sep=";")

    # Load dictionary of ligand xyz coordinates
    ligand_xyzs_path = os.path.join(input_dir, "xyz", "ligands_xyzs.xyz")
    ligand_xyzs = load_ligand_xyz(ligand_xyzs_path)

    # Get denticity-specific Dataframe from tmQMg-L
    if denticity == "monodentate":
        df = get_monodentate(df_fingerprints, df_misc)
        output_csv_path = os.path.join(args.output_dir, "tmQMg-L_mono.csv")
        output_txt_path = os.path.join(args.output_dir, "tmQMg-L_mono.txt")
        core_element = "Li"
    elif denticity == "bidentate":
        df = get_bidentate(df_fingerprints, df_misc)
        output_csv_path = os.path.join(args.output_dir, "tmQMg-L_bi.csv")
        output_txt_path = os.path.join(args.output_dir, "tmQMg-L_bi.txt")
        core_element = "Ir"

    # Sample for debugging
    if args.debug:
        df = df[0:10].copy(deep=True)

    start = time.time()
    with mp.Pool(mp.cpu_count()) as pool:
        mols = pool.starmap(
            ligand_xyz_to_mol, [(row, df_stable) for row in df.iterrows()]
        )
    end = time.time()
    print(f"Total time: {end-start}")
    print(mols[0:10])

    mol_objects = [x[0] if x is not None else None for x in mols]
    connect_id = [x[1] if x is not None else None for x in mols]
    canon_smile = [Chem.MolToSmiles(x) if x is not None else None for x in mol_objects]
    df["custom_mol"] = mol_objects
    df["connect_id"] = connect_id
    df["Canonical SMILES"] = canon_smile
    df.dropna(inplace=True)
    df["enriched_mol"] = df.apply(
        attach_dummy_atom_to_coordinating_atoms,
        axis=1,
        element=core_element,
        joint=True,
    )
    df.dropna(subset=["enriched_mol"], inplace=True)
    df["enriched_mol"] = df["enriched_mol"].swifter.apply(prune_num_atoms, num=4)
    df.dropna(subset=["enriched_mol"], inplace=True)
    mols = df["enriched_mol"]
    mols = [Chem.RemoveHs(x) for x in mols]
    [Chem.RemoveStereochemistry(x) for x in mols]
    [x.RemoveAllConformers() for x in mols]
    smi = [Chem.MolToSmiles(x) for x in mols]
    df["enriched_smiles"] = smi

    with open(os.path.join(__location__, "discarded_ligand_ids.txt"), "r") as f:
        discard_ligands = f.read().splitlines()

    # Drop all rows that represent manually discarded ligands.
    df = df.drop(df[df["name"].isin(discard_ligands)].index)
    df = df.rename(
        columns={
            "name": "tmQMg-L ligand ID",
            "Canonical SMILES": "Canonical SMILES",
            "connect_id": "Connection IDs",
            "enriched_smiles": "Enriched SMILES",
        }
    )
    df = df[
        ["tmQMg-L ligand ID", "Canonical SMILES", "Connection IDs", "Enriched SMILES"]
    ]

    df.to_csv(output_csv_path, index=False)
    # Created enriched SMILES training file
    df["Enriched SMILES"].to_csv(output_txt_path, header=False, index=False)

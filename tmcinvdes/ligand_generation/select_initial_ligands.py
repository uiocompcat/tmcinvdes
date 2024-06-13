import argparse

# from itertools import repeat
# from functools import partial
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
    get_id,
    get_monodentate,
    get_stable_occ,
    load_ligand_xyz,
    prune_num_atoms,
)


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
    connect_ids = get_id(row, df_stable)
    #
    for elem in connect_ids:
        connect_id = elem[0]
        stable_occ = get_stable_occ(row["name"], df_stable)

        xyz = ligand_xyzs[stable_occ]

        if isinstance(row["charge"], np.integer):
            charge = row["charge"].item()
        else:
            charge = row["charge"]

        # Check if we can get smiles from tmQMg-L SMILES
        try:
            m = Chem.MolFromSmiles(row["smiles"])
            if not m:
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

        # Check for carbene, that should not be a radical.
        elif m.GetSubstructMatch(
            Chem.MolFromSmarts("[#6&v2H0,#6&v3H0]")
        ) or m.GetSubstructMatch(Chem.MolFromSmarts("[#14&v2H0]")):
            print("babel")
            mol = x2m_tm.get_mol_pure_babel(xyz, charge=charge)
        elif any(x in row["smiles"] for x in ["As", "Se", "p", "P"]):
            mol = x2m_tm.get_mol_pure_babel(xyz, charge=charge)
        else:
            mol = x2m_tm.get_mol(xyz, charge=charge)
        if not mol:
            return None

        # If there are still radicals present,
        for atom in mol.GetAtoms():
            if atom.GetNumRadicalElectrons() == 1:
                # filtered_stuff['radicals_after_mol'].append(row['name'])
                return None

        connect_atom = mol.GetAtomWithIdx(connect_id)
        id = connect_atom.GetIdx()

        valence = mol.GetAtomWithIdx(id).GetExplicitValence()

        if connect_atom.GetSymbol() == "P":
            if valence > 4:
                return None

        if connect_atom.GetSymbol() == "C":
            if valence > 3:
                return None
            if connect_atom.GetFormalCharge() == -1:
                return None

    return (mol, connect_ids, stable_occ)


if __name__ == "__main__":
    # global ligand_xyzs
    # ligand_xyzs = load_ligand_xyz()

    args = parse_args()
    denticity = args.denticity

    args.output_dir.mkdir(exist_ok=True, parents=True)

    input_dir = os.path.abspath(args.input_dir)
    df_desc = pd.read_csv(os.path.join(input_dir, "ligands_descriptors.csv"), sep=";")
    df_fingerprints = pd.read_csv(
        os.path.join(input_dir, "ligands_fingerprints.csv"), sep=";"
    )
    df_misc = pd.read_csv(os.path.join(input_dir, "ligands_misc_info.csv"), sep=";")
    df_stable = pd.read_csv(os.path.join(input_dir, "stable.csv"), sep=";")
    ligand_xyzs_path = os.path.join(input_dir, "xyz", "ligands_xyzs.xyz")
    ligand_xyzs = load_ligand_xyz(ligand_xyzs_path)

    # Get denticity-specific Dataframe from tmQMg-L
    if denticity == "monodentate":
        df = get_monodentate(df_fingerprints, df_misc)
    elif denticity == "bidentate":
        df = get_bidentate(df_fingerprints, df_misc)

    # Sample for debugging
    if args.debug:
        df = df[0:10].copy(deep=True)

    start = time.time()
    with mp.Pool(mp.cpu_count()) as pool:
        # df_stable = pd.read_csv(os.path.join(input_dir, "stable.csv"), sep=";")
        # partial_ligand_xyz_to_mol = partial(ligand_xyz_to_mol, df_stable=df_stable)
        mols = pool.starmap(
            ligand_xyz_to_mol, [(row, df_stable) for row in df.iterrows()]
        )
    end = time.time()
    print(f"Total time: {end-start}")
    print(mols[0:10])

    mol = [x[0] if x is not None else None for x in mols]
    connect_id = [x[1] if x is not None else None for x in mols]
    df["custom_mol"] = mol
    df["connect_id"] = connect_id
    df.dropna(inplace=True)
    df["enriched_mol"] = df.apply(
        attach_dummy_atom_to_coordinating_atoms,
        axis=1,
        element="Li",
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

    # These ligands were manually discarded:
    manual_discard = [
        "ligand5249-0",
        "ligand9550-0",
        "ligand11200-0",
        "ligand11246-0",
        "ligand13615-0",
        "ligand17783-0",
        "ligand19660-0",
        "ligand24383-0",
    ]

    # Maybe additional discards?
    reproduced_excess = [
        "ligand1994-1",
        "ligand3740-0",
        "ligand4236-0",
        "ligand13618-0",
        "ligand20893-0",
        "ligand30481-0",
        "ligand30860-0",
        "ligand31507-0",
        "ligand31574-0",
    ]
    reproduced_deficit = ["ligand25836-0"]  # Failed to reproduce this one.
    #

    # Drop all rows that represent manually discarded ligands.
    df = df.drop(df[df["name"].isin(manual_discard)].index)
    df = df.drop(df[df["name"].isin(reproduced_excess)].index)

    df.to_csv(args.output_dir / f"train_{denticity}.csv", index=False)

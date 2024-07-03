"""Process and screen unconditionally generated ligands to extract only novel
ligands.

The script also screens out generated ligands that are invalid.
The training set is necessary as a reference to determine which ligands are novel relative to the
training set.

Usage:

python screen_generated_ligands.py -d denticity -i input_dir -o output_dir -t training_set_dir -x test

python screen_generated_ligands.py -d bidentate \
         -i ../datasets/03_uncond-raw50k-generated_ligands \
         -o ../datasets/04_uncond-novel \
         -t ../datasets/01_tmQMg-L-training_sets \
         -x test
"""

import argparse
import os
import re
import sys
from pathlib import Path
from typing import Callable

import pandas as pd
from rdkit import Chem

from tmcinvdes.ligand_generation.utils import (
    compare_dataframes,
    get_neighbors_bidentate,
    get_smiles_donor_id,
    process_substitute_attachment_points,
    process_substitute_attachment_points_bidentate,
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
        default="bidentate",
        help="Select one of the two denticity modes supported.",
    )
    parser.add_argument(
        "--input_dir",
        "-i",
        type=Path,
        required=True,
        help="""Input directory with the raw, unprocessed files of generated ligands as enriched
        SMILES.""",
    )
    parser.add_argument(
        "--output_dir",
        "-o",
        type=Path,
        required=True,
        help="Output directory for the initial screen of the generated ligands.",
    )
    parser.add_argument(
        "--training_set_dir",
        "-t",
        type=Path,
        required=True,
        help="""Directory of training sets used to train the model used for generating the ligands to
        be screened.""",
    )
    parser.add_argument(
        "--xtent",
        "-x",
        choices=["full", "test", "demo"],
        default="test",
        help="""The extent to which the process should run.
         - In the "full" case, the process runs in its entirety and writes the output results to
           file.
         - In the "test" case, the process runs in its entirety and compares the results to the
           output already existing in the designated output file.
         - In the "demo" case, the process only runs to a limited extent (n ~ 100 instances) and
           if applicable compares the limited results to the corresponding instances already
           existing in the designated output file.""",
    )
    return parser.parse_args(arg_list)


def count_atoms_hfull(mol: Chem.rdchem.Mol) -> int:
    """Count the number of atoms in the molecule, including hydrogens.

    Args:
        mol (Chem.rdchem.Mol): molecule whose number of atoms is being counted.

    Returns:
        atom_count (int): number of atoms in molecule, including hydrogen atoms.
    """
    rwmol = Chem.RWMol(mol)
    atom_count = 0
    molwh = Chem.AddHs(rwmol)
    for _ in molwh.GetAtoms():
        atom_count += 1
    assert atom_count == molwh.GetNumAtoms()
    return atom_count


def count_atoms_hless(mol: Chem.rdchem.Mol) -> int:
    """Count the number of heavy atoms in the molecule.

    Args:
        mol (Chem.rdchem.Mol): molecule whose number of atoms is being counted.

    Returns:
        atom_count (int): number of atoms in molecule, excluding hydrogen atoms.
    """
    rwmol = Chem.RWMol(mol)
    atom_count = 0
    molwo = Chem.RemoveHs(rwmol)
    for _ in molwo.GetAtoms():
        atom_count += 1
    assert atom_count == molwo.GetNumAtoms()
    return atom_count


# Simple regex functions to count substring patterns in enriched SMILES as strings:


def count_berylium(smile):
    regex = r"Be"
    matches = re.findall(regex, smile)
    return len(matches)


def count_lithium(smile):
    regex = r"Li"
    matches = re.findall(regex, smile)
    return len(matches)


def count_iridium(smile):
    regex = r"Ir"
    matches = re.findall(regex, smile)
    return len(matches)


def count_left_li(smile):
    regex = r"\[Li\]<-"
    matches = re.findall(regex, smile)
    return len(matches)


def count_left_be(smile):
    regex = r"\[Be\]="
    matches = re.findall(regex, smile)
    return len(matches)


def count_left_carbene(smile):
    regex = r"\[Be\]=\[C\]"
    matches = re.findall(regex, smile)
    return len(matches)


def count_left_sylene(smile):
    regex = r"\[Be\]=\[Si\]"
    matches = re.findall(regex, smile)
    return len(matches)


def count_dative(smile):
    regex1 = r"->"
    matches1 = re.findall(regex1, smile)
    regex2 = r"<-"
    matches2 = re.findall(regex2, smile)
    return len(matches1) + len(matches2)


def count_dative_bidentate_ir(smile):
    regex = r"->\[Ir\]<-"
    matches = re.findall(regex, smile)
    return len(matches)


def count_right_dative(smile):
    regex = r"->"
    matches = re.findall(regex, smile)
    return len(matches)


def get_coordination_env(mol: Chem.rdchem.Mol, connection_ids: list) -> list:
    """Get the coordination environments of the connection atoms.

    Args:
        mol (Chem.rdchem.Mol): the Mol object for neighborhood reference.
        connection_ids (list): the IDs of the connection atoms of the ligand.

    Returns:
        list: the string(s) corresponding to the coordination environment of each connection atom.
    """
    envs = []
    for connection_id in connection_ids:
        atom = mol.GetAtomWithIdx(connection_id)
        neighbors = atom.GetNeighbors()
        elements = [
            f"([{a.GetSymbol()}])"
            for a in neighbors
            if (
                a.GetSymbol() != "Li"
                and a.GetSymbol() != "Be"
                and a.GetSymbol() != "Ir"
            )
        ]

        # Distinguish carbenes
        for neigh in neighbors:
            if neigh.GetSymbol() == "Be":
                carb = "_carb"
                break
        else:
            carb = ""

        elements.insert(0, f"{atom.GetSymbol()}{carb}")

        pattern = "".join(elements)
        envs.append(pattern)
    return envs


def process_enriched_smiles_monodentate(
    enriched_smilesx: list,
    atom_counting_mode: str,
    count_atoms_fun: Callable,
    process_substitutions_fun: Callable,
) -> dict:
    """Process enriched SMILES strings into Mol objects and canonical SMILES,
    or else discard them. Monodentate version.

    Args:
        enriched_smilesx (list): list of the validated subset of enriched SMILES, Mol objects and
        connection atoms.
        atom_counting_mode (str): the mode of counting atoms per ligand.
        count_atoms_fun (Callable): the function used to count atoms.
        process_substitutions_fun (Callable): the function used to process the substitutions from
        the enriched SMILES encoding.

    Returns:
        dict: a map of the successfully processed ligands, with canonical SMILES as key for easy
        comparison between generated and training ligands.
    """
    processed_dict = {}

    for smile, mol, connection_id in enriched_smilesx:
        atom = mol.GetAtomWithIdx(connection_id)
        neighbors = atom.GetNeighbors()

        if atom.GetSymbol() == "P":
            val = atom.GetTotalValence()
            if val > 4:
                continue
        flag = False
        for neigh in neighbors:
            val = neigh.GetTotalValence()
            if val > 5:
                flag = True
        if flag:
            continue
        envs = get_coordination_env(mol, [connection_id])

        try:
            new_mol, new_connection_ids = process_substitutions_fun(mol)
            if isinstance(
                new_connection_ids, int
            ):  # Always for monodentates, left for later refactor
                new_connection_ids = [new_connection_ids]
            canon_smiles, hless_mapped_ids = get_smiles_donor_id(new_mol)
            try:
                canon_connection_ids = [hless_mapped_ids[x] for x in new_connection_ids]
            except Exception:
                if connection_id == 1:
                    canon_connection_ids = [connection_id - 1]
                else:
                    print(smile)
                    print(connection_id)
                    # Maybe right sided enrichment?
                    print(
                        f"Maybe try connection ID as: {len([a for a in mol.GetAtoms()])-1}"
                    )
                    sys.exit()
            atom_count = count_atoms_fun(mol)
            processed_dict[canon_smiles] = {
                "mol": new_mol,
                "Enriched SMILES": smile,
                "Connection IDs": canon_connection_ids,
                atom_counting_mode: atom_count,
                "Coordination environment": envs,
            }
        except Exception:
            continue
    return processed_dict


def process_enriched_smiles_bidentate(
    enriched_smilesx: list,
    atom_counting_mode: str,
    count_atoms_fun: Callable,
    process_substitutions_fun: Callable,
) -> dict:
    """Process enriched SMILES strings into Mol objects and canonical SMILES,
    or else discard them.

    Args:
        enriched_smilesx (list): list of enriched SMILES strings.
        atom_counting_mode (str): the mode of counting atoms per ligand.
        count_atoms_fun (Callable): the function used to count atoms.
        process_substitutions_fun (Callable): the function used to process the substitutions from
        the enriched SMILES encoding.

    Returns:
        dict: a map of the successfully processed ligands, with canonical SMILES as key for easy
        comparison between generated and training ligands.
    """
    processed_dict = {}
    for smile in enriched_smilesx:
        mol = Chem.MolFromSmiles(smile)
        connection_ids = get_neighbors_bidentate(mol)
        if connection_ids is None:
            continue
        else:
            envs = get_coordination_env(mol, connection_ids)
            envs = sorted(envs)
        try:
            mol, connection_ids = process_substitutions_fun(mol)
            if isinstance(connection_ids, int):
                connection_ids = [connection_ids]
            canon_smiles, hless_mapped_ids = get_smiles_donor_id(mol)
            try:
                canon_connection_ids = [hless_mapped_ids[x] for x in connection_ids]
            except Exception:
                continue
            atom_count = count_atoms_fun(mol)
            # Mark for deletion ligands where the connection IDs are not unambiguous:
            if canon_smiles in processed_dict:
                if (
                    processed_dict[canon_smiles]["Connection IDs"]
                    != canon_connection_ids
                ):
                    processed_dict[canon_smiles] = {
                        "mol": None,
                        "Enriched SMILES": smile,
                        "Connection IDs": "delete me",
                    }
                else:
                    continue
            else:
                processed_dict[canon_smiles] = {
                    "mol": mol,
                    "Enriched SMILES": smile,
                    "Connection IDs": canon_connection_ids,
                    atom_counting_mode: atom_count,
                    "Coordination environment": envs,
                }
        except Exception:
            continue
    return processed_dict


def remove_ambiguously_connected_ligands(processed_dict: dict) -> dict:
    """Remove processed ligands where multiple connection modes were present,
    creating ambiguity about which atoms would coordinate with the metal core
    of a TMC.

    Args:
        processed_dict (dict): map of successfully process ligands, including ligands marked for
        deletion.

    Returns:
        dict: map of successfully process ligands, after removing ligands marked for deletion.
    """
    delkeys = []
    for smiles in processed_dict.keys():
        if processed_dict[smiles]["Connection IDs"] == "delete me":
            delkeys.append(smiles)
    for key in delkeys:
        del processed_dict[key]
    return processed_dict


def get_unique_in_order(xs: list) -> list:
    """Remove non-unique instances of an object after its first occurrence.

    Args:
        xs (list): list of objects with potential duplicates.

    Returns:
        list: only the first occurrence of each unique object.
    """
    series = pd.Series(xs)
    unique_ordered = series.unique().tolist()
    return unique_ordered


def validate_enriched_smiles_monodentate(enriched_smilesx: list) -> list:
    """Validate enriched SMILES that represent monodentate ligands.

    This removes enriched SMILES which fall outside the expected form of generated ligands, at a
    string level of testing.

    Args:
        enriched_smilesx (list): the enriched SMILES to validate.

    Returns:
        list: the validated subset of enriched SMILES, Mol objects and connection atoms.
    """
    enriched_smilesx.sort(key=lambda x: Chem.MolFromSmiles(x).GetNumAtoms())
    smiles_mols = [
        tuple([x, Chem.MolFromSmiles(x)])
        for x in enriched_smilesx
        if (Chem.MolFromSmiles(x))
    ]
    be_smart = Chem.MolFromSmarts("[Be]")
    li_smart = Chem.MolFromSmarts("[Li]")
    connection_atoms = []
    valid_mols = []
    valid_smilesx = []
    for smiles, mol in smiles_mols:
        id = mol.GetSubstructMatch(be_smart)
        id2 = mol.GetSubstructMatch(li_smart)
        if id or id2:
            atom_idx = mol.GetAtomWithIdx(0)
            neighbor_atom = atom_idx.GetNeighbors()
            if not neighbor_atom:
                continue
            neighbor_idx = neighbor_atom[0].GetIdx()
            connection_atoms.append(neighbor_idx)
            valid_mols.append(mol)
            valid_smilesx.append(smiles)
        # else: # Just drop
        #    connection_atoms.append(None)
    return zip(valid_smilesx, valid_mols, connection_atoms)


def validate_enriched_smiles_bidentate(enriched_smilesx: list) -> list:
    """Validate enriched SMILES that represent bidentate ligands.

    This removes enriched SMILES which fall outside the expected form of generated ligands, at a
    string level of testing.

    Args:
        enriched_smilesx (list): the enriched SMILES to validate.

    Returns:
        list: the validated subset of enriched SMILES.
    """
    validated_enriched_smilesx = []
    for smile in enriched_smilesx:
        count_bid = count_dative_bidentate_ir(smile)
        count_dat = count_dative(smile)
        count_ir = count_iridium(smile)
        if count_ir == 1 and count_dat == 2 and count_bid == 1:
            validated_enriched_smilesx.append(smile)
    return validated_enriched_smilesx


def screen_generated_ligands(
    df_input: pd.DataFrame,
    df_train: pd.DataFrame,
    denticity: str,
    atom_counting_mode: str,
    validate_enriched_smiles_fun: Callable,
    count_atoms_fun: Callable,
    process_substitutions_fun: Callable,
) -> pd.DataFrame:
    """Screen generated ligands for validity and novelty relative to the
    training set.

    Args:
        df_input (pd.DataFrame): the dataframe with the unprocessed generated ligand set as
        enriched SMILES.
        df_train (pd.DataFrame): the dataframe with generating model's training set of ligands as
        enriched SMILES.
        atom_counting_mode (str): the mode of counting atoms, depending on denticity.

    Returns:
        pd.DataFrame: the output dataframe with the final set of screened and processed ligands.
    """
    cols = [
        atom_counting_mode,
        "Canonical SMILES",
        "Connection IDs",
        "Enriched SMILES",
        "Coordination environment",
    ]
    unique_generated = get_unique_in_order(df_input[0].values.tolist())

    unique_train = get_unique_in_order(df_train[0].values.tolist())
    # Remove generated enriched SMILES that are verbatim string matches of training set instances:
    unique_generated = [x for x in unique_generated if x not in unique_train]
    validated_generated = validate_enriched_smiles_fun(unique_generated)
    # No need to validate the training set of Enriched SMILES.
    if denticity == "monodentate":  # Shared code for monodentates diverges
        processed_generated = process_enriched_smiles_monodentate(
            validated_generated,
            atom_counting_mode,
            count_atoms_fun,
            process_substitutions_fun,
        )
    elif denticity == "bidentate":
        processed_generated = process_enriched_smiles_bidentate(
            validated_generated,
            atom_counting_mode,
            count_atoms_fun,
            process_substitutions_fun,
        )
        processed_generated = remove_ambiguously_connected_ligands(processed_generated)
        processed_train = process_enriched_smiles_bidentate(
            unique_train, atom_counting_mode, count_atoms_fun, process_substitutions_fun
        )
        # No need to remove ambiguously connected ligands in the training set.
        # Remove generated enriched SMILES that are structural matches of training set instances:
        for train_canon_smiles in processed_train.keys():
            if train_canon_smiles in processed_generated:
                del processed_generated[train_canon_smiles]
        processed_generated_sorted = {
            k: v
            for k, v in sorted(
                processed_generated.items(), key=lambda x: x[1][atom_counting_mode]
            )
        }
        processed_generated = processed_generated_sorted

    data_dict = {col: [] for col in cols}
    for smiles, vals in processed_generated.items():
        data_dict[atom_counting_mode].append(vals[atom_counting_mode])
        data_dict["Canonical SMILES"].append(smiles)
        data_dict["Connection IDs"].append(sorted(vals["Connection IDs"]))
        data_dict["Enriched SMILES"].append(vals["Enriched SMILES"])
        data_dict["Coordination environment"].append(vals["Coordination environment"])
    df_output = pd.DataFrame(data_dict, columns=cols)
    if atom_counting_mode == "Heavy atom count":
        df_output = df_output[df_output["Heavy atom count"] > 3]
        df_output = df_output.reset_index(drop=True)
    return df_output


if __name__ == "__main__":
    args = parse_args()
    denticity = args.denticity
    xtent = args.xtent
    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output_dir)
    train_set_dir = os.path.abspath(args.training_set_dir)

    # Modify parameters based on denticity:
    if denticity == "monodentate":
        input_file = os.path.join(input_dir, "uncond_mono-raw50k.txt")
        output_file = os.path.join(output_dir, "uncond_mono-raw50k-novel.csv")
        train_set_file = os.path.join(train_set_dir, "tmQMg-L_mono.txt")
        atom_counting_mode = "Heavy atom count"
        validate_enriched_smiles_fun = validate_enriched_smiles_monodentate
        count_atoms_fun = count_atoms_hless
        process_substitutions_fun = process_substitute_attachment_points
    elif denticity == "bidentate":
        input_file = os.path.join(input_dir, "uncond_bi-raw50k.txt")
        output_file = os.path.join(output_dir, "uncond_bi-raw50k-novel.csv")
        train_set_file = os.path.join(train_set_dir, "tmQMg-L_bi.txt")
        atom_counting_mode = "Atom count"
        validate_enriched_smiles_fun = validate_enriched_smiles_bidentate
        count_atoms_fun = count_atoms_hfull
        process_substitutions_fun = process_substitute_attachment_points_bidentate

    # These files should have no column header, each is just a single column of enriched SMILES.
    df_input = pd.read_csv(input_file, header=None)
    df_train = pd.read_csv(train_set_file, header=None)

    # Modify process based on xtent:
    if xtent == "demo":  # Constrain how many instances to process.
        pass  # Will probably drop this mode, especially for non-order-preserving processes.
    elif xtent == "test":
        df_expect = pd.read_csv(output_file)
        df_expect["Connection IDs"] = [
            eval(x) for x in df_expect["Connection IDs"].values.tolist()
        ]
        df_output = screen_generated_ligands(
            df_input,
            df_train,
            denticity,
            atom_counting_mode,
            validate_enriched_smiles_fun,
            count_atoms_fun,
            process_substitutions_fun,
        )
        perfect_match = df_output.equals(df_expect)
        row_accuracy = compare_dataframes(df_output, df_expect, atom_counting_mode)
        if perfect_match or row_accuracy == 1.0:
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
        df_output = screen_generated_ligands(
            df_input,
            df_train,
            denticity,
            atom_counting_mode,
            validate_enriched_smiles_fun,
            count_atoms_fun,
            process_substitutions_fun,
        )
        df_output.to_csv(output_file, index=False, header=True)

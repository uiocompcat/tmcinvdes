"""Provide utility functions related to handling ligand data."""

import ast
import json
import os
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem


def make_dict_xyz(ligand_xyzs_path: Path):
    """Create dict of from the single xyz ligand file.

    The keys in the dict will be the name of the TMC for which the
    ligand xyz is found.
    """
    xyzs = {}
    with open(ligand_xyzs_path, "r") as fh:
        for xyz in fh.read().split("\n\n"):
            xyzs[xyz.split("\n")[1]] = xyz

    with open("ligands_dict_xyz.json", "w") as f:
        json.dump(xyzs, f)
    print("Succesfully created the dict of stable ligand complexes")


def load_ligand_xyz(ligand_xyzs_path: Path):
    if not os.path.isfile("ligands_dict_xyz.json"):
        print("Dict with ligand xyz coordinates does not excist. Creating it now.")
        make_dict_xyz(ligand_xyzs_path)

    # load ligand xyz dict
    with open("ligands_dict_xyz.json", "r") as f:
        xyzs = json.load(f)

    return xyzs


def get_monodentate(df, df2, charge=0):
    "Get df of neutral monodentates"

    df2["charge"] = df["charge"]
    monodentate = df2[(df["n_metal_bound"] == 1) & (df["charge"] == charge)]
    return monodentate


def get_bidentate(df, df2, charge=0):
    "Get dataframe of neutral monodentates and bidentates."
    df2["charge"] = df["charge"]
    bi_mask = (df["n_metal_bound"] == 2) & (df["n_dentic_bound"] == 2)
    charge_mask = df2["charge"] == charge
    monodentate = df2[(bi_mask) & charge_mask]
    return monodentate


def get_stable_occ(name: str, df_stable: pd.DataFrame):
    "Get string of most stable occurence label for a ligand"
    stable_oc = df_stable[df_stable["name"] == name]["stable_occurrence_name"].item()
    return stable_oc


def get_connection_ids(row, df_stable):
    "Extract the connection id of ligand as list"
    stable_oc = get_stable_occ(row["name"], df_stable)
    res = row["metal_bond_node_idx_groups"]
    ocs = ast.literal_eval(res)
    connection_ids = ocs[stable_oc]
    return connection_ids


def _smarts_filter(mol: Chem.rdchem.Mol, connect_ids: list[int]):
    """Helper function to find various exceptions to the connection points.

    Currently carbenes and sylenes are detected here. For these we use a
    Be= connection point.
    """
    type_matches = []

    # Check for carbene
    match_c = mol.GetSubstructMatches(Chem.MolFromSmarts("[#6&v2H0,#6&v3H0]"))
    # Check for sylene
    match_si = mol.GetSubstructMatches(Chem.MolFromSmarts("[Si]"))

    for id in connect_ids:
        if any(id in match for match in match_c):
            type_matches.append("carbene")
            continue

        if any(id in match for match in match_si):
            if mol.GetAtomWithIdx(id).GetNumRadicalElectrons() == 2:
                type_matches.append("sylene")
                continue

        # If none of the matches were found, we set generic label.
        type_matches.append("generic")

    return type_matches


def attach_dummy_atom_to_coordinating_atoms(row, element="Ir", joint=False):
    """Function to attach dummy substitutes to the ligands to create dataset
    for JT-VAE."""
    mol = row["custom_mol"]
    connect_id = row["connect_id"]

    # Extract the connection IDs.
    if len(connect_id) > 1:
        ids = [int(x[0]) for x in connect_id]
    else:
        ids = [int(connect_id[0][0])]

    # To ensure radicals are present in the Mol object, we need to sanitize.
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        print("Sanitation failed with error:")
        print(e)
        return None

    # Start editing Mol.
    emol = Chem.RWMol(mol)
    emol.BeginBatchEdit()

    type_matches = _smarts_filter(mol, ids)

    mapper = {
        "carbene": {"bond_type": Chem.BondType.DOUBLE, "substitute_element": "Be"},
        "sylene": {"bond_type": Chem.BondType.DOUBLE, "substitute_element": "Be"},
        "generic": {"bond_type": Chem.BondType.DATIVE, "substitute_element": element},
    }

    alternate_idxs = []
    idx_generic = False
    first_match = True
    for connect_id, type_match in zip(ids, type_matches):
        if type_match == "generic":
            if joint:
                if first_match:
                    idx_generic = emol.AddAtom(
                        Chem.Atom(mapper[type_match]["substitute_element"])
                    )
                    first_match = False
            else:
                idx_generic = emol.AddAtom(
                    Chem.Atom(mapper[type_match]["substitute_element"])
                )
            emol.AddBond(int(connect_id), idx_generic, mapper[type_match]["bond_type"])
        else:
            idx = emol.AddAtom(Chem.Atom(mapper[type_match]["substitute_element"]))
            alternate_idxs.append(idx)
            emol.AddBond(int(connect_id), idx, mapper[type_match]["bond_type"])

    # Decide whether to connect everything on the same Iridium.
    if joint:
        if alternate_idxs and len(type_matches) > 1:
            if idx_generic:
                for alt_id in alternate_idxs:
                    emol.AddBond(alt_id, idx_generic, Chem.BondType.DATIVE)
            else:
                idx_generic = emol.AddAtom(
                    Chem.Atom(mapper["generic"]["substitute_element"])
                )
                for alt_id in alternate_idxs:
                    emol.AddBond(alt_id, idx_generic, Chem.BondType.DATIVE)

    emol.CommitBatchEdit()
    mol = emol.GetMol()
    try:
        Chem.SanitizeMol(mol)
        return mol
    except Exception:
        return None


def get_smiles_donor_id(mol: Chem.rdchem.Mol):
    """When the Mol object is written to file, the atom ordering can change.

    This function gets the SMILES and maps the ordering of the original
    Mol object atoms to match the new ordering in the SMILES.
    """
    mol = Chem.RemoveHs(mol)
    smi = Chem.MolToSmiles(mol)
    order = eval(mol.GetProp("_smilesAtomOutputOrder"))
    mapped_ids = {}
    for i in range(0, len(mol.GetAtoms())):
        mapped_id = np.where(np.array(order) == i)[0][0]
        mapped_ids[i] = mapped_id
    return smi, mapped_ids


def prune_num_atoms(mol: Chem.rdchem.Mol, num: int = None):
    if isinstance(mol, str):
        try:
            mol = Chem.MolFromSmiles(mol)
            if not mol:
                return None
        except Exception as e:
            print("Mol from SMILES failed with error: ")
            print(e)
            return None
    try:
        num_atoms = Chem.RemoveHs(mol).GetNumAtoms()
    except Exception as e:
        print("Removing hydrogens failed with error:")
        print(e)
        return None

    # Remove the weird SIs.
    smart = Chem.MolFromSmarts("[Li]<-[Si]")
    if mol.GetSubstructMatch(smart):
        return None

    if num_atoms > num:
        try:
            Chem.SanitizeMol(mol)
            Chem.Kekulize(mol)
            smi = Chem.MolToSmiles(mol)
            if "." in smi:
                return None

            return mol
        except Exception:
            return None
    else:
        return None


def single_atom_remover(mol: Chem.rdchem.Mol, idx: int) -> Chem.rdchem.Mol:
    """Function that removes an atom at specified idx.

    Args:
        mol (Chem.rdchem.Mol): The Mol from which to remove atom.
        idx (int): ID of atom to remove from the input Mol.

    Returns:
        The Mol with the atom removed.
    """
    res = Chem.RWMol(mol)
    res.BeginBatchEdit()
    res.RemoveAtom(idx)
    res.CommitBatchEdit()
    Chem.SanitizeMol(res)
    return res.GetMol()


def process_substitute_attachment_points(mol: Chem.rdchem.Mol) -> tuple:
    """Remove substitute attachment points. This decodes a monodentate ligand,
    so that the Mol can output a decoded SMILES.

    Args:
        mol (Chem.rdchem.Mol): the encoded Mol of a monodentate ligand.

    Returns:
        Decoded Mol object and its connection ID.
    """
    # Determine which attachment.
    substitute_smarts = Chem.MolFromSmarts("[Be,Li]")

    matches = mol.GetSubstructMatches(substitute_smarts)
    # If there are several matches for monodentates, then molecule should be discarded.
    connect_id = None
    if len(matches) > 1:
        new_mol = None
    elif not matches:
        new_mol = None
    else:
        match_atom = mol.GetAtomWithIdx(matches[0][0])

        # We remove the substitute here. Since Li or Be is always at 0, the coordinating atom
        # will get id 0 after removal.
        try:
            new_mol = single_atom_remover(mol, matches[0][0])
            connect_id = 0
        except Exception as e:
            print("Single atom remover failed with error: ")
            print(e)
            return None, None

        # If Be, then the carbon Be was connected to should be a Carbene.
        if match_atom.GetSymbol() == "Be":
            # 2 radical atoms to create the carbene.
            new_mol.GetAtomWithIdx(0).SetNumRadicalElectrons(2)

            # Ensure that there are no hydrogens on the Carbon atom.
            new_mol.GetAtomWithIdx(0).SetNoImplicit(True)
            new_mol.GetAtomWithIdx(0).SetNumExplicitHs(0)

    return new_mol, connect_id


def process_substitute_attachment_points_bidentate(mol: Chem.rdchem.Mol) -> tuple:
    """Remove substitute attachment points. This decodes a bidentate ligand, so
    that the Mol can output a decoded SMILES.

    Args:
        mol (Chem.rdchem.Mol): the encoded Mol of a bidentate ligand.

    Returns:
        Decoded Mol object and its connection IDs.
    """
    # Determine the Ir substitute atom idx.
    substitute_smarts = Chem.MolFromSmarts("[Ir]")
    matches = mol.GetSubstructMatches(substitute_smarts)

    # If there are several matches for Ir, then we need to discard the mol.
    connection_ids = None
    if len(matches) > 1:
        new_mol = None
    elif not matches:
        new_mol = None
    else:
        # Get the neighbors of Ir.
        neighbors = mol.GetAtomWithIdx(matches[0][0]).GetNeighbors()

        # For bidentates, there should be exactly 2 neighbors.
        if len(neighbors) != 2:
            print("There are not two neighbors to the Ir")
            return None, None
        try:
            tmp_mol = single_atom_remover(mol, matches[0][0])
        except Exception as e:
            print("Single atom remover failed with error: ")
            print(e)
            return None, None

        # Loop through neighbors to create connection IDs list.
        connection_ids = []
        for n in neighbors:
            id = n.GetIdx()
            # If Be, we get the ID later, after Be removal.
            if n.GetSymbol() == "Be":
                continue

            # Since we are removing an atom, any atoms with idx higher than the idx of the Ir atom,
            # will decrease by one.
            # We therefore need to check for this to set the right connection ID.
            if id > matches[0][0]:
                new_id = id - 1
                connection_ids.append(new_id)
            else:
                connection_ids.append(id)
        # Finally, if one of the neighbors to the Ir is a Be, we need to remove this as well.
        be_match = tmp_mol.GetSubstructMatches(Chem.MolFromSmarts("[Be]"))

        if be_match:
            res = Chem.RWMol(tmp_mol)
            res.BeginBatchEdit()

            for match in be_match:
                carbene_neighbor = res.GetAtomWithIdx(match[0]).GetNeighbors()[0]
                carbene_neighbor_idx = carbene_neighbor.GetIdx()

                # We need to explicitely ensure that the carbon atom now is a carbene with 2
                # radical electrons and 0 hydrogens.
                carbene_neighbor.SetNumRadicalElectrons(2)
                carbene_neighbor.SetNoImplicit(True)
                carbene_neighbor.SetNumExplicitHs(0)

                res.RemoveAtom(match[0])

                # The idx of the new carbene will decrease with number of Be removed, if the
                # idx of the carbene is larger than the Be idx.
                idx_decrease = sum(i[0] < carbene_neighbor_idx for i in be_match)
                if connection_ids and idx_decrease == 0 and len(be_match) == 1:
                    connection_ids[0] -= 1
                carbene_neighbor_idx = carbene_neighbor_idx - idx_decrease
                connection_ids.append(carbene_neighbor_idx)
            res.CommitBatchEdit()
            Chem.SanitizeMol(res)
            tmp_mol = res.GetMol()
        new_mol = tmp_mol
    return new_mol, connection_ids


def get_neighbors_bidentate(mol: Chem.rdchem.Mol) -> list:
    """For the Mol object of a bidentate ligand constructed from an encoded
    SMILES string, find the IDs of the connection atoms before removing the
    extra atoms of the encoding.

    Args:
        mol (Chem.rdchem.Mol): from encoded SMILES string.

    Returns:
        list: the IDs of the connection atoms in the Mol object.
    """
    connection_ids = (
        []
    )  # Atom IDs of the connection atoms before removing the encoding atom(s).
    substitute_smarts = Chem.MolFromSmarts("[Ir]")
    matches = mol.GetSubstructMatches(substitute_smarts)
    if len(matches) == 1:
        atom = mol.GetAtomWithIdx(matches[0][0])
        neighbors = atom.GetNeighbors()
        assert len(neighbors) > 0
        for neighbor in neighbors:
            if neighbor.GetSymbol() not in ["Be", "Li", "Ir"]:
                connection_ids.append(neighbor.GetIdx())
            else:
                nextneighbors = neighbor.GetNeighbors()
                assert len(neighbors) > 0
                for nextneighbor in nextneighbors:
                    if nextneighbor.GetSymbol() not in ["Be", "Li", "Ir"]:
                        connection_ids.append(nextneighbor.GetIdx())
                        break
        return connection_ids
    else:
        return None


def read_file(file_name: Path, num_mols: int) -> list:
    """Read SMILES from file and return as list."""
    mols = []
    with open(file_name, "r") as file:
        for i, smiles in enumerate(file):
            mols.append(smiles)
            if i == num_mols:
                break
    return mols


def compare_dataframes(
    df_output: pd.DataFrame,
    df_expect: pd.DataFrame,
    atom_counting_mode: str = "Heavy atom count",
) -> float:
    """Compare the reproduced and expected dataframe and calculate accuracy in
    terms of identical rows.

    Args:
        df_output (pd.DataFrame): the reproduced dataframe that would be the output written to file
        if not testing.
        df_expect (pd.DataFrame): the target output read back into a dataframe from file.
        atom_counting_mode (str): the mode of counting atoms per ligand.

    Returns:
        float: accuracy as a number between 0.0 and 1.0 representing the proportion of overlapping
        identical rows between the dataframes being compared.
    """
    df_output["Connection IDs"] = df_output["Connection IDs"].apply(str)
    df_expect["Connection IDs"] = df_expect["Connection IDs"].apply(str)
    df_output["Coordination environment"] = df_output["Coordination environment"].apply(
        str
    )
    df_expect["Coordination environment"] = df_expect["Coordination environment"].apply(
        str
    )
    rows_intersect = pd.merge(
        df_output,
        df_expect,
        how="inner",
        on=[
            atom_counting_mode,
            "Decoded SMILES",
            "Connection IDs",
            "Encoded SMILES",
            "Coordination environment",
        ],
    )
    rows_union = pd.merge(
        df_output,
        df_expect,
        how="outer",
        on=[
            atom_counting_mode,
            "Decoded SMILES",
            "Connection IDs",
            "Encoded SMILES",
            "Coordination environment",
        ],
    )
    if rows_intersect.equals(rows_union):
        return 1.0
    row_accuracy = float(rows_intersect.shape[0]) / rows_union.shape[0]
    return row_accuracy


def get_elem_counter_hless(mol):
    elem_counter = defaultdict(int)
    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        if sym != "H":
            elem_counter[sym] += 1
    return dict(sorted(elem_counter.items(), key=lambda item: item[0]))


def check_structural_match(mol1: Chem.rdchem.Mol, mol2: Chem.rdchem.Mol) -> bool:
    """Attempt in various way to match two Mol objects.

    Args:
        mol1 (Chem.rdchem.Mol): the first Mol to match.
        mol2 (Chem.rdchem.Mol): the second Mol to match.

    Returns:
        True if matched, otherwise False.
    """
    assert mol1 is not None and mol2 is not None
    # Compare initial input Mols.
    if mol1.HasSubstructMatch(mol2) and mol2.HasSubstructMatch(mol1):
        return True
    # Try to ensure both Mols are stripped of implicit Hydrogens
    try:
        mol1_hless = Chem.RemoveHs(mol1)
        mol2_hless = Chem.RemoveHs(mol2)
        if mol1_hless.HasSubstructMatch(mol2_hless) and mol2_hless.HasSubstructMatch(
            mol1_hless
        ):
            return True
    except Exception:
        pass
    try:
        mol1_hless = Chem.RemoveAllHs(mol1)
        mol2_hless = Chem.RemoveAllHs(mol2)
        if mol1_hless.HasSubstructMatch(mol2_hless) and mol2_hless.HasSubstructMatch(
            mol1_hless
        ):
            return True
    except Exception:
        pass
    try:
        mol1_hfull = Chem.AddHs(mol1_hless)
        mol2_hfull = Chem.AddHs(mol2_hless)
        if mol1_hfull.HasSubstructMatch(mol2_hfull) and mol2_hfull.HasSubstructMatch(
            mol1_hfull
        ):
            return True
    except Exception:
        pass
    return False


def compare_encoded_smiles_structures(
    df_output: pd.DataFrame, df_expect: pd.DataFrame, denticity: str = "monodentate"
) -> float:
    """Compare two sets of encoded SMILES of ligands to find their proportional
    overlap, i.e., accuracy, not as verbatim strings but as Mol objects.

    Args:
        df_output (pd.DataFrame): encoded SMILES, the reproduced output from the process step.
        df_expect (pd.DataFrame): encoded SMILES, the expected output read from file.
        denticity (str, optional): the denticity of all the ligands.

    Returns:
        The accuracy score based on matched Mol structures.
    """
    output_mol_dict = defaultdict(list)
    expect_mol_dict = defaultdict(list)
    if denticity == "monodentate":
        process_function = process_substitute_attachment_points
    elif denticity == "bidentate":
        process_function = process_substitute_attachment_points_bidentate

    for i, row in df_expect.iterrows():
        encoded_smiles = row["Encoded SMILES"]
        mol = Chem.MolFromSmiles(encoded_smiles)
        new_mol, _ = process_function(mol)
        if new_mol is None:
            continue
        query_elems = str(get_elem_counter_hless(new_mol))
        expect_mol_dict[query_elems].append(
            {"Encoded SMILES": encoded_smiles, "mol": new_mol}
        )

    for i, row in df_output.iterrows():
        ligand_id = row["tmQMg-L ligand ID"]
        encoded_smiles = row["Encoded SMILES"]
        mol = Chem.MolFromSmiles(encoded_smiles)
        new_mol, _ = process_function(mol)
        if new_mol is None:
            continue
        query_elems = str(get_elem_counter_hless(new_mol))
        output_mol_dict[query_elems].append(
            {
                "Encoded SMILES": encoded_smiles,
                "mol": new_mol,
                "tmQMg-L ligand ID": ligand_id,
            }
        )

    output_to_expect_strucmatch_dict = {}
    matched_orig_encoded_smiles = set()

    for element_key, reprod_vals in output_mol_dict.items():
        orig_vals = expect_mol_dict[element_key]
        for reprod_val in reprod_vals:
            reprod_encoded_smiles = reprod_val["Encoded SMILES"]
            reprod_mol = reprod_val["mol"]
            ligand_id = reprod_val["tmQMg-L ligand ID"]
            for orig_val in orig_vals:
                orig_encoded_smiles = orig_val["Encoded SMILES"]
                if orig_encoded_smiles in matched_orig_encoded_smiles:
                    continue
                orig_mol = orig_val["mol"]
                if check_structural_match(reprod_mol, orig_mol):
                    output_to_expect_strucmatch_dict[reprod_encoded_smiles] = {
                        "Original encoded SMILES": orig_encoded_smiles,
                        "tmQMg-L ligand ID": ligand_id,
                    }
                    matched_orig_encoded_smiles.add(orig_encoded_smiles)
                    break

    unmatched_orig_encoded_smiles = set()
    for i, row in df_expect.iterrows():
        encoded_smiles = row["Encoded SMILES"]
        if encoded_smiles not in matched_orig_encoded_smiles:
            unmatched_orig_encoded_smiles.add(encoded_smiles)

    unmatched_reprod_encoded_smiles = set()
    for i, row in df_expect.iterrows():
        encoded_smiles = row["Encoded SMILES"]
        if encoded_smiles not in output_to_expect_strucmatch_dict.keys():
            unmatched_reprod_encoded_smiles.add(encoded_smiles)
    structure_accuracy = len(matched_orig_encoded_smiles) / (
        len(matched_orig_encoded_smiles)
        + len(unmatched_orig_encoded_smiles)
        + len(unmatched_reprod_encoded_smiles)
    )

    return structure_accuracy

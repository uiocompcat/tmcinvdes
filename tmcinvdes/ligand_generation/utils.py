import ast
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from openbabel import openbabel
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


def get_smiles_openbabel_hannes(xyz: str):
    """Gets the SMILES string of a given xyz structure using Hannes' Open Babel
    method that was used to obtain smiles in tmQMg-L. NB! This method will give
    you uncharged SMILES for charged ligands. It also cannot place atomic
    charged like xyz2mol.

    Arguments:
        xyz (str): The xyz structure.

    Returns:
        str: The SMILES string.
    """
    # setup converter for reading xyz
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("xyz", "mdl")

    # setup molecule
    mol = openbabel.OBMol()
    obConversion.ReadString(mol, xyz)

    # set converter for getting smiles
    obConversion.SetInAndOutFormats("mdl", "smi")

    return obConversion.WriteString(mol).strip()


def get_monodentate(df, df2, charge=0):
    "Get df of neutral monodentates"

    df2["charge"] = df["charge"]
    monodentate = df2[(df["n_metal_bound"] == 1) & (df["charge"] == charge)]
    return monodentate


def get_bidentate(df, df2, charge=0):
    "Get dataframe of neutral monodentates and bidentates."
    df2["charge"] = df["charge"]
    # mono_mask = (df["n_metal_bound"] == 1) & (df["n_dentic_bound"] == 1)
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


def _smarts_filter(mol, connect_ids):
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

    # Extract the binding ids
    if len(connect_id) > 1:
        ids = [int(x[0]) for x in connect_id]
    else:
        ids = [int(connect_id[0][0])]

    # To ensure radicals are present in the mol object, we need to sanitize
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        print("Sanitation failed with error:")
        print(e)
        return None

    # Start editing mol
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

    # Decide wether to connect everything on the same Irridium
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
        # print(Chem.MolToSmiles(mol))
        return mol
    except Exception:
        return None


def get_smiles_donor_id(mol):
    """When mol object is written to file the atom ordering can change.

    This function gets the SMILES and maps the ordering of the original
    mol object atoms to match the new ordering in the SMILES.
    """
    mol = Chem.RemoveHs(mol)
    smi = Chem.MolToSmiles(mol)
    order = eval(mol.GetProp("_smilesAtomOutputOrder"))
    mapped_ids = {}
    for i in range(0, len(mol.GetAtoms())):
        mapped_id = np.where(np.array(order) == i)[0][0]
        mapped_ids[i] = mapped_id
    return smi, mapped_ids


def prune_num_atoms(mol, num=None):
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

    # Remove the weird SIs
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


def single_atom_remover(mol, idx):
    """Function that removes an atom at specified idx.

    Args:
        mol (Chem.rdchem.Mol): The Mol to remove substruct on
        idx (int): idx of atom to remove from the input Mol

    Returns:
        Chem.rdchem.Mol: The ouput Mol with the atom removed
    """
    res = Chem.RWMol(mol)
    res.BeginBatchEdit()
    res.RemoveAtom(idx)
    res.CommitBatchEdit()
    Chem.SanitizeMol(res)
    return res.GetMol()


def process_substitute_attachment_points(mol):
    "Remove substitute attachment points"

    # Determine which attachment. Not working for Ir yet.
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

        # We remove the substitute here. Since Li or Be is always 0, the coordinating atom
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

            # Ensure that there are no hydrogens on the Carbon atom
            new_mol.GetAtomWithIdx(0).SetNoImplicit(True)
            new_mol.GetAtomWithIdx(0).SetNumExplicitHs(0)

    return new_mol, connect_id


def process_substitute_attachment_points_bidentate(mol):
    # Determine the Ir substitute atom idx
    substitute_smarts = Chem.MolFromSmarts("[Ir]")
    matches = mol.GetSubstructMatches(
        substitute_smarts
    )  # TL: What is the structure of matches and its elements and subelements?

    # If there are several matches for Ir, then we need to discard the mol
    connection_ids = None
    if len(matches) > 1:
        new_mol = None
    elif not matches:
        new_mol = None
    else:
        # Get Ir - Atom object
        # match_atom = mol.GetAtomWithIdx()
        # Get neighbors to Ir
        neighbors = mol.GetAtomWithIdx(matches[0][0]).GetNeighbors()

        # For bidentates, there should be exactly 2 neighbors
        if len(neighbors) != 2:
            print("There are not two neighbors to the Ir")
            return None, None
        try:
            tmp_mol = single_atom_remover(mol, matches[0][0])
        except Exception as e:
            print("Single atom remover failed with error: ")
            print(e)
            return None, None

        # Loop through neighbors to create connection_id list
        connection_ids = []
        for n in neighbors:
            id = n.GetIdx()
            # If Be, we get the id later after Be removal.
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
                # idx of the carbene i larger than the Be idx. This is a bit of a hack.
                idx_decrease = sum(i[0] < carbene_neighbor_idx for i in be_match)
                # TL: { hack to fix the hack{
                if connection_ids and idx_decrease == 0 and len(be_match) == 1:
                    connection_ids[0] -= 1  # always one?
                    # print("Decreasing index of non-Be neighbour by 1")
                # }
                carbene_neighbor_idx = carbene_neighbor_idx - idx_decrease
                connection_ids.append(carbene_neighbor_idx)
            res.CommitBatchEdit()
            Chem.SanitizeMol(res)
            tmp_mol = res.GetMol()
        new_mol = tmp_mol
    return new_mol, connection_ids


def get_neighbors_bidentate(mol: Chem.rdchem.Mol) -> list:
    connection_ids = (
        []
    )  # atom IDs of the connection atoms before removing the enrichment atom(s).
    substitute_smarts = Chem.MolFromSmarts("[Ir]")
    matches = mol.GetSubstructMatches(substitute_smarts)
    if len(matches) == 1:
        atom = mol.GetAtomWithIdx(matches[0][0])
        neighbors = atom.GetNeighbors()
        assert len(neighbors) > 0
        for neighbor in neighbors:
            connection_ids.append(neighbor.GetIdx())
        return connection_ids
    else:
        return None


def read_file(file_name, num_mols):
    """Read smiles from file and return Mol generator."""
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
            "Canonical SMILES",
            "Connection IDs",
            "Enriched SMILES",
            "Coordination environment",
        ],
    )
    rows_union = pd.merge(
        df_output,
        df_expect,
        how="outer",
        on=[
            atom_counting_mode,
            "Canonical SMILES",
            "Connection IDs",
            "Enriched SMILES",
            "Coordination environment",
        ],
    )
    if rows_intersect.equals(rows_union):
        return 1.0
    row_accuracy = float(rows_intersect.shape[0]) / rows_union.shape[0]
    return row_accuracy

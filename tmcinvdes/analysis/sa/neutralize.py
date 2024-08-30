"""Functionality to neutralize molecules."""

import copy
import json
from pathlib import Path
from typing import List, Tuple

from rdkit import Chem
from rdkit.Chem.rdmolops import FastFindRings

from .sascorer import calculateScore

_neutralize_reactions = None

package_directory = Path(__file__).parent.resolve()


def read_neutralizers(name="neutralize") -> List[Tuple[Chem.Mol, Chem.Mol]]:
    filename = package_directory / f"{name}.json"
    with open(filename) as json_file:
        reactions = json.load(json_file)
        neutralize_reactions = []
        for reaction in reactions:
            reaction["name"]
            r_s = reaction["reactant"]
            p_s = reaction["product"]
            r = Chem.MolFromSmarts(r_s)
            p = Chem.MolFromSmiles(p_s, False)
            assert r is not None and p is not None, "Either R or P is None"
            neutralize_reactions.append((r, p))
    return neutralize_reactions


def neutralize_smiles(smiles: List[str]) -> List[str]:
    """Neutralize a set of SMILES.

    :param smiles: a list of SMILES
    """
    assert type(smiles) == list

    charged_molecules = [Chem.MolFromSmiles(s) for s in smiles]
    neutral_molecules = neutralize_molecules(charged_molecules)
    return [Chem.MolToSmiles(m) for m in neutral_molecules]


def neutralize_molecules(population):
    """Neutralize a set of molecules.

    :param list charged_molecules: list of (possibly) charged molecules
    :return: list of neutral molecules
    """
    assert type(population.molecules) == list
    global _neutralize_reactions
    if _neutralize_reactions is None:
        _neutralize_reactions = read_neutralizers()

    for individual in population.molecules:
        mol = copy.deepcopy(individual.rdkit_mol_sa)
        mol.UpdatePropertyCache()
        Chem.rdmolops.FastFindRings(mol)
        assert mol is not None
        for reactant_mol, product_mol in _neutralize_reactions:
            while mol.HasSubstructMatch(reactant_mol):
                rms = Chem.ReplaceSubstructs(mol, reactant_mol, product_mol)
                if rms[0] is not None:
                    mol = rms[0]
        # https://github.com/rdkit/rdkit/issues/2216 to get valid ring information and implicit valences
        mol.UpdatePropertyCache()
        FastFindRings(mol)
        individual.neutral_rdkit_mol = mol


if __name__ == "__main__":
    s_q = "c1ccccc1C(C(=O)[O-])c2ccccc2"
    s_n = neutralize_smiles([s_q])

    print("Q: ", s_q, calculateScore(Chem.MolFromSmiles(s_q)))
    print("N: ", s_n, calculateScore(Chem.MolFromSmiles(s_n[0])))

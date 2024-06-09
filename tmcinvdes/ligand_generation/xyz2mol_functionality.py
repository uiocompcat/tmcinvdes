"""Module doing xyz -> mol conversions using either build in RDKIt
fucntionality or OpenBabel."""

import os
import subprocess
import sys
import time
import uuid

from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

# NB this is for using a local installation of xyz2mol.
# This is now build into RDKit and is not really needed.
module_path = os.path.abspath("/home/magstr/git/xyz2mol/")
if module_path not in sys.path:
    sys.path.append(module_path)
import xyz2mol as x2m

TRANSITION_METALS = [
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "La",
    "Ni",
    "Cu",
    "Zn",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
]


def read_file(file_name, num_mols):
    """Read smiles from file."""
    mols = []
    with open(file_name, "r") as file:
        for i, smiles in enumerate(file):
            mols.append(smiles.rstrip())
            if i == num_mols:
                break
    return mols


def setChiralTagAndOrder(atom, chiralTag, chiralPermutation=None):
    """Sets the chiral tag of an atom and the permutation order of attached
    ligands. These tags are used in RDKit embedding into 3D structures.

    Args:
    atom (Chem.Atom): Atom for which to set the chiral tag/permutation order properties
    chiralTag (Chem.rdchem.ChiralType, optional): Chiral Tag of Metal Atom. Defaults to None.
    permutationOrder (int, optional): Permutation order of ligands. Defaults to None.
    """
    atom.SetChiralTag(chiralTag)
    if chiralPermutation:
        atom.SetIntProp("_chiralPermutation", chiralPermutation)


def shell(cmd, shell=False):
    "Utility function to run shell commands"
    if shell:
        p = subprocess.Popen(
            cmd,
            shell=True,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    else:
        cmd = cmd.split()
        p = subprocess.Popen(
            cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )

    output, err = p.communicate()
    return output


def get_connect_mol(xyz):
    "Get adjacency matrix"
    raw_mol = Chem.MolFromXYZBlock(xyz)
    mol = Chem.Mol(raw_mol)
    rdDetermineBonds.DetermineConnectivity(mol, useHueckel=True)

    return mol


def get_mol(xyz, charge):
    """Get mol object from xyz coordinates.

    openbabel can give better connectivity for mols where xyz2mol gives
    fragments hense the if statement.
    """
    # Get simply mol with connectivity, bond orders not determined.
    mol = get_connect_mol(xyz)

    # If theres fragment in smiles, try to use babel to get connectivity.
    if "." in Chem.MolToSmiles(mol, canonical=False):
        mol = get_mol_babel(xyz, charge)
        return mol
    else:
        # This function fails for some elements. Therefore a try statement is needed.
        try:
            rdDetermineBonds.DetermineBonds(mol, charge=charge, useHueckel=True)
            return mol
        except Exception as e:
            print(f"Failed xyz: {e}")
            return None


def get_mol_pure_babel(xyz, charge):
    "Use commandline version of openbabel to get final mol"

    # Create unique string. In this was the functions can be called in parallel and write in
    # the same directory without intertwining
    unique_identifier = uuid.uuid4()

    raw_mol = Chem.MolFromXYZBlock(xyz)
    Chem.MolToXYZFile(raw_mol, f"{unique_identifier}_test.xyz")
    cmd = (
        "obabel -ixyz "
        + f" {unique_identifier}_test.xyz -O {unique_identifier}_test.sdf"
    )
    _ = shell(cmd, shell=False)
    time.sleep(0.1)  # time needed to finish writing file
    suppl = Chem.SDMolSupplier(
        f"{unique_identifier}_test.sdf", removeHs=False, sanitize=False
    )
    # Remove written files
    shell(f"rm {unique_identifier}_test.sdf", shell=False)
    shell(f"rm {unique_identifier}_test.xyz", shell=False)

    # openbabel can provide resonance structures. Not sure how to enable this.
    # I always only see one mol object in suppl
    mols = [x for x in suppl]
    mol = mols[0]

    # We need to check if there are still fragments. If we could not
    # find mol without fragments return None.
    if mol and "." not in Chem.MolToSmiles(mol):
        return mol
    else:
        print("Failed pure babel")
        return None

    return mol


def get_mol_babel(xyz, charge):
    """Get connectivity from openbabel and then use xyz2mol to get the full mol
    objet."""

    raw_mol = Chem.MolFromXYZBlock(xyz)

    unique_identifier = uuid.uuid4()
    Chem.MolToXYZFile(raw_mol, f"{unique_identifier}_test.xyz")
    cmd = (
        "obabel -ixyz "
        + f"{unique_identifier}_test.xyz -O {unique_identifier}_test.sdf"
    )
    _ = shell(cmd, shell=False)
    time.sleep(0.1)  # time needed to finish writing file
    suppl = Chem.SDMolSupplier(
        f"{unique_identifier}_test.sdf", removeHs=False, sanitize=False
    )
    shell(f"rm {unique_identifier}_test.sdf", shell=False)
    shell(f"rm {unique_identifier}_test.xyz", shell=False)
    mols = [x for x in suppl]
    mol = mols[0]

    charged_fragments = True
    quick = True

    atoms = [a.GetAtomicNum() for a in mol.GetAtoms()]
    adjacent_matrix = Chem.GetAdjacencyMatrix(mol)

    # Define new molecule template from atoms
    new_mol = x2m.get_proto_mol(atoms)

    try:
        # reconstruct the molecule from adjacent matrix, atoms and total charge
        new_mols = x2m.AC2mol(
            new_mol, adjacent_matrix, atoms, charge, charged_fragments, quick
        )
    except Exception as e:
        print(f"AC2mol failed, e= {e}")
        return None
    if new_mols and "." not in Chem.MolToSmiles(new_mols[0]):
        new_mol = new_mols[0]
        # Optional embedding of molecule. Will slow down function significantly
        # rdDistGeom.EmbedMolecule(new_mol, useRandomCoords=True, maxAttempts=10)
        return new_mol
    else:
        print("Failed babel")
        return None

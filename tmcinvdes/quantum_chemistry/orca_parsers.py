"""Parse ORCA output/logfile.

Usage:

python orca_parsers.py -i orca.out
"""

import argparse
from pathlib import Path
from pprint import pprint
from typing import List

from .cm5 import get_cm5


def get_args(arg_list=None):
    """Parse CLI arguments.

    Args:
        arg_list: Automatically obtained from the commandline if provided.
        If no arguments are given but default arguments are defined, the latter are used.

    Returns:
        parser.parse_args(arg_list)(Namespace): Dictionary-like class that contain the driver arguments.
    """
    parser = argparse.ArgumentParser(description="Parse ORCA output file/log.")
    parser.add_argument(
        "--input_file",
        "-i",
        type=Path,
        required=True,
        help="Path to input file, the log of the ORCA job. Typically `orca.out`.",
    )
    return parser.parse_args(arg_list)


def normal_termination(lines: List[str]) -> bool:
    """Check if ORCA terminated normally."""
    for line in reversed(lines):
        if line.strip() == "****ORCA TERMINATED NORMALLY****":
            return True
    return False


def read_final_sp_energy(lines: List[str]) -> float:
    """Read final single point energy from ORCA output."""
    for line in reversed(lines):
        if "FINAL SINGLE POINT ENERGY" in line:
            return float(line.split()[-1])
    return None


def read_homo_lumo(lines) -> float:
    """Read the homo lumo gap."""
    start_spin_down = None
    start_spin_up = None
    start_spin = None
    previous_line = [""]
    for i, line in enumerate(reversed(lines)):
        if "SPIN UP ORBITALS" in line:
            start_spin_up = i - 1
            break
        if "SPIN DOWN ORBITALS" in line:
            start_spin_down = i - 1
        if "ORBITAL ENERGIES" in line and not (start_spin_up or start_spin_down):
            start_spin = i - 3
            break

    if start_spin_up:
        for line in lines[-start_spin_up:]:
            split_line = line.split()
            occ = int(float(split_line[1]))
            if occ != 1 and occ != 2:
                up_unocc_energy = float(split_line[2])
                up_occ_energy = float(previous_line[2])
                break
            previous_line = split_line

        for line in lines[-start_spin_down:]:
            split_line = line.split()
            occ = int(float(split_line[1]))
            if occ != 1 and occ != 2:
                down_unocc_energy = float(split_line[2])
                down_occ_energy = float(previous_line[2])
                break
            previous_line = split_line
        highest_unoccupied = min([down_unocc_energy, up_unocc_energy])
        lowest_occupied = max([down_occ_energy, up_occ_energy])
    else:
        for line in lines[-start_spin:]:
            split_line = line.split()
            occ = int(float(split_line[1]))
            if occ != 1 and occ != 2:
                highest_unoccupied = float(split_line[2])
                lowest_occupied = float(previous_line[2])
                break
            previous_line = split_line

    homo_lumo_energy = highest_unoccupied - lowest_occupied

    return homo_lumo_energy


def read_opt_structure(lines: List[str]) -> tuple:
    """Read optimized structure from ORCA output.

    Args:
        lines (List[str]): ORCA output as a list of strings.

    Returns:
        tuple: (atoms, coords)
    """
    rev_start_idx = 2  # Just a dummy value
    for i, line in enumerate(reversed(lines)):
        if "CARTESIAN COORDINATES (ANGSTROEM)" in line:
            rev_start_idx = i - 1
            break

    atoms = []
    coords = []
    for line in lines[-rev_start_idx:]:
        line = line.strip()
        if len(line) > 0:
            atom, x, y, z = line.split()
            atoms.append(atom)
            coords.append([float(x), float(y), float(z)])
        else:
            break
    return atoms, coords


def repackage_xyz(name: str, atoms: list, coords: list) -> str:
    """Re-package the opt_structure output into an XYZ block as a string.

    Args:
        name (str): name of the TMC/structure
        atoms (list): list of atomic symbols.
        coords (list): list of XYZ coordinates.

    Returns:
        An XYZ block as a string.
    """
    lines = []
    lines.append(str(len(atoms)))
    lines.append(name)
    for atom, coord in zip(atoms, coords):
        terms = [atom]
        for c in coord:
            terms.append(str(c))
        assert len(terms) == 4
        lines.append(" ".join(terms))
    return "\n".join(lines)


def get_orca_results(
    lines: List[str], properties: List[str] = ["electronic_energy"]
) -> dict:
    """Read results from ORCA output.

    Args:
        lines (List[str]): ORCA output as a list of strings.
        properties (List[str], optional): Properties. Defaults to ["electronic_energy"].


    Returns:
        dict: Results dictionary.
    """
    assert isinstance(lines, list), "Input lines must be a list of strings"

    results = {}
    reader = {
        "electronic_energy": read_final_sp_energy,
        "opt_structure": read_opt_structure,
        "homo_lumo_energy": read_homo_lumo,
        "cm5": get_cm5,
    }

    for property in properties:
        results[property] = reader[property](lines)

    return results


def main():
    args = get_args()
    with open(args.input_file, "r") as f:
        lines = f.readlines()
    res = get_orca_results(
        lines,
        properties=["electronic_energy", "opt_structure", "homo_lumo_energy", "cm5"],
    )
    pprint(res)


if __name__ == "__main__":
    main()

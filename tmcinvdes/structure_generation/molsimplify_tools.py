"""Shared tools for use in structure generation scripts, to standardize
processes and criteria."""

import os
import shlex
import shutil
import subprocess

import uxtbpy


def mybash(cmd):
    """Function to run a bash command, waiting 8 min to complete."""
    cmd = shlex.split(cmd)
    try:
        cp = subprocess.run(
            cmd, shell=False, capture_output=True, text=True, timeout=480
        )
        stout = str(cp.stdout)
    except subprocess.TimeoutExpired:
        stout = "Process timed out at 8 minutes. Attempting to kill process."
        print(stout)
    return stout


def get_existing_ligand_names():
    """Return the ligand names present in the molSimplify library.

    Returns:
        list: The list of existing ligand names.
    """
    out = str(
        subprocess.run(
            ["molsimplify", *["-h", "liganddict"]],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ).stdout
    )

    return [_.strip() for _ in out.replace("}", "").split("\\n")[8:-1]]


def add_ligand_from_xyz(xyz: str, name: str, connection_ids: list):
    """Add a ligand to the molSimplify library.

    Args:
        xyz (str): The XYZ structure of the ligand.
        name (str): The name of the ligand to be used as identifier in the molSimplify library.
        connection_ids (list): The list of connection indicies (0-indexed).
    """
    if name in get_existing_ligand_names():
        print(
            'Cannot add ligand because the name "'
            + name
            + '" already exists! Please choose a different name.'
        )
        return

    # write molecule to temporary file
    with open("temp_mol.xyz", "w") as fh:
        fh.write(xyz)

    parameters = [
        "-ligadd temp_mol.xyz",
        "-ligname " + name,
        "-ligcon " + ",".join([str(_ + 1) for _ in connection_ids]),
        "-skipANN True",
    ]
    # subprocess.run(
    #     ["molsimplify", *parameters], stdout=subprocess.PIPE, stderr=subprocess.PIPE
    # )
    mybash(" ".join(["molsimplify", *parameters]))

    # clean up temporary file
    os.remove("temp_mol.xyz")


def add_ligand_from_smiles(
    smiles: str, name: str, connection_ids: list, pre_checked: bool = False
):
    """Add a ligand to the molSimplify library.

    Args:
        smiles (str): The (decoded) SMILES string of the ligand.
        name (str): The name of the ligand to be used as identifier in the molSimplify library.
        connection_ids (list): The list of connection indicies (0-indexed).
        pre_checked (bool): Whether already checked for presence in pre-existing local ligands.
    """
    if not pre_checked:
        if name in get_existing_ligand_names():
            print(
                'Cannot add ligand because the name "'
                + name
                + '" already exists! Please choose a different name.'
            )
            return

    parameters = [
        "-ligadd " + smiles,
        "-ligname " + name,
        "-ligcon " + ",".join([str(_ + 1) for _ in connection_ids]),
        "-skipANN True",
    ]

    mybash(" ".join(["molsimplify", *parameters]))


def build_tmc(geometry: str, metal_center: str, ligand_names: list):
    """Assemble a transition metal complex (TMC) with ligands.

    Args:
        geometry (str): The metal complex geometry (li, tpl, sqp, thd, spy, tbp, oct, tpr, pbp, sqap).
        metal_center (str): The metal center element identifier, i.e., atomic symbol.
        ligand_names (list): List of the ligand names.

    Returns:
        str: The XYZ structure of the built metal complex.
    """
    # Remove molSimplify operation directory in case files from previous run are present.
    shutil.rmtree("Runs/", ignore_errors=True)

    parameters = [
        "-skipANN True",
        "-core " + metal_center,
        "-geometry " + geometry,
        "-lig " + ",".join(ligand_names),
        "-ligocc " + ",".join(["1" for _ in ligand_names]),
        "-name run",
        "-keepHs yes",
    ]
    mybash(" ".join(["molsimplify", *parameters]))

    with open("Runs/run/run/run.xyz") as fh:
        xyz = fh.read()
        return xyz


def calculate_cis_and_trans_energy(
    metal_center: str, ligand_names: list, charge: int
) -> dict:
    """Calculates the xtb energies for the cis and trans coordination
    geometries for a bidentate square planar complex.

    Args:
        metal_center (str): The metal center element identifier, i.e., atomic symbol.
        ligand_names (list): The list of ligand names as they are stored in the
            molSimplify library. Expects ligand names without the '_flipped' suffix.
        charge (int): The overall charge of the TMC.

    Returns:
        The energies of cis and trans isomerisms of the otherwise equivalent homoleptic TMCs.
    """
    # set up xtb runner

    xtb_runner = uxtbpy.XtbRunner(output_format="dict")

    # Assume every occurrence of a ligand is spelled out:
    assert len(ligand_names) == 2

    trans_xyz = build_tmc("sqp", metal_center, ligand_names)
    ligand_names_flipped = ligand_names
    ligand_names_flipped[1] += "_flipped"
    cis_xyz = build_tmc("sqp", metal_center, ligand_names_flipped)

    cis_energy = xtb_runner.run_xtb_from_xyz(
        cis_xyz, parameters=["--opt tight", "--norestart", f"-c {charge}"]
    )["energy"]

    trans_energy = xtb_runner.run_xtb_from_xyz(
        trans_xyz, parameters=["--opt tight", "--norestart", f"-c {charge}"]
    )["energy"]

    return {"cis_energy": cis_energy, "trans_energy": trans_energy}


def calculate_fac_and_mer_energy(metal_center: str, ligand_names: list, charge: int):
    """Calculates the xtb energies for the fac and mer coordination geometries
    for a bidentate octahedral complex.

    Args:
        metal_center (str): The metal center element identifier, i.e., atomic symbol.
        ligand_names (list): The list of ligand names as they are stored in the
            molSimplify library. Expects ligand names without the '_flipped' suffix.
        charge (int): The overall charge of the TMC.

    Returns:
        dict: Result dict containing the FAC and MER energies.
    """
    # Set up xtb runner.

    xtb_runner = uxtbpy.XtbRunner(output_format="dict")

    assert len(ligand_names) == 3

    mer_xyz = build_tmc("oct", metal_center, ligand_names)
    ligand_names_flipped = ligand_names
    ligand_names_flipped[1] += "_flipped"
    fac_xyz = build_tmc("oct", metal_center, ligand_names_flipped)

    fac_energy = xtb_runner.run_xtb_from_xyz(
        fac_xyz, parameters=["--opt tight", "--norestart", f"-c {charge}"]
    )["energy"]

    mer_energy = xtb_runner.run_xtb_from_xyz(
        mer_xyz, parameters=["--opt tight", "--norestart", f"-c {charge}"]
    )["energy"]

    return {"fac_energy": fac_energy, "mer_energy": mer_energy}


if __name__ == "__main__":
    # Generic example of usage:
    add_ligand_from_smiles("[O]CCC[O]", "test", [0, 4])
    print(build_tmc("oct", "Fe", ["test", "test", "test"]))

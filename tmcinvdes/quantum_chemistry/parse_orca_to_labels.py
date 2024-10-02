"""Parse ORCA logs to create .csv with labeled files.

ORCA calculations have been performed on TMCs whose initial structures were generated using
molSimplify.

This script takes as input a directory of successful ORCA calculations' logfiles and returns the
parsed, labeled results as a .csv-file.

The filenames of each ORCA logfile for each unique homoleptic TMC are assumed to reflect the
ligand ID, and in the case of bidentate TMCs, also the isomerism.

Usage: (from root directory of repository)

python -m tmcinvdes.quantum_chemistry.parse_orca_to_labels -d bidentate \
                                    -i tmp_results/orca_out-uncond_bi-min10k-TMC/ \
                                    -o datasets/07_uncond-labeled/uncond_bi-min10k-labeled.csv \
                                    -s datasets/05_uncond-minXk/uncond_bi-min10k.csv \
                                    -x test
"""

import argparse
import os
from pathlib import Path

import morfeus
import pandas as pd
from parse import parse
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from scipy import constants

from tmcinvdes.ligand_generation.utils import (
    process_substitute_attachment_points,
    process_substitute_attachment_points_bidentate,
)
from tmcinvdes.quantum_chemistry.orca_parsers import get_orca_results, repackage_xyz
from tmcinvdes.utils import compare_dataframes

Eh_to_eV = constants.physical_constants["Hartree energy in eV"][0]


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
        default="tmp_results/orca_out-uncond_bi-min10k-TMC/",
        required=True,
        help="Input directory with only successful ORCA .out files, each named after the TMC ID",
    )
    parser.add_argument(
        "--output_file",
        "-o",
        type=Path,
        default="datasets/07_uncond-labeled/uncond_bi-min10k-labeled.csv",
        required=True,
        help="Output filename for .csv, same prefix used for .xyz files.",
    )
    parser.add_argument(
        "--reference_optimization",
        "-r",
        type=Path,
        default="datasets/11_cond-optimized/cond_mono-sampled_optimized.csv",
        required=False,
        help="""If ligands were optimized, provide an additional reference input file of optimized
        ligands used to construct input TMCs. Conditionally required argument.""",
    )
    parser.add_argument(
        "--smiles_reference",
        "-s",
        type=Path,
        default="datasets/05_uncond-minXk/uncond_bi-min10k.csv",
        required=False,
        help="""If input TMCs were constructed from unconditionally generated bidenate ligands,
        provide an additional reference input file to the decoded SMILES. Conditionally required
        argument.""",
    )
    parser.add_argument(
        "--type",
        "-t",
        choices=["unoptimized", "optimized"],
        default="unoptimized",
        help="""The type of ligands used to build TMCs that were labeled with ORCA. If unoptimized,
        there is less prior information to include.""",
    )
    parser.add_argument(
        "--xtent",
        "-x",
        choices=["full", "test"],
        default="test",
        help="""The extent to which the process should run.
         - In the "full" case, the process runs in its entirety and writes the output results to
           file.
         - In the "test" case, the process runs in its entirety and compares the results to the
           output already existing in the designated output file.""",
    )
    return parser.parse_args(arg_list)


def get_g_parameter(xyz: str, num_ligands: int) -> float:
    """Obtain the solid angle G parameter of a single ligand from a homoleptic
    TMC geometry.

    Args:
        xyz (str): the geometry of a homoleptic TMC.
        num_ligands (int): the number of identical ligands in the TMC.

    Returns:
        The G parameter.
    """
    xyz_lines = [line.rstrip() for line in xyz.rstrip("\n\n").split("\n")]
    num_atoms = int(xyz_lines[0])
    assert xyz_lines[2][:2] == "Ir"
    assert len(xyz_lines) == 2 + num_atoms  #:

    offset = (num_atoms - 1) / num_ligands
    assert 2 * (int(offset)) == int(2 * offset)
    offset = int(offset)

    xyz_lines = xyz_lines[: 3 + offset]
    new_num_atoms = 1 + offset
    xyz_lines[0] = str(new_num_atoms)
    xyz = "\n".join(xyz_lines)

    # read xyz into MORFEUS
    tmp_filename = "temp.xyz"
    with open(tmp_filename, "w") as fh:
        fh.write(xyz)
    elements, coordinates = morfeus.read_xyz(tmp_filename)
    solid_angle_calculator = morfeus.SolidAngle(elements, coordinates, 1)
    g_parameter = solid_angle_calculator.G
    return g_parameter


def add_ligand_descriptors(
    df: pd.DataFrame, denticity, optimized: bool = False, df_smiles: pd.DataFrame = None
) -> pd.DataFrame:
    """Add ligand descriptors of individual ligands.

    Specifically, this function evaluates and adds the descriptors log P and solid angle
    G parameter to the output dataframe.

    Args:
        df (pd.DataFrame): dataframe to be extended.
        denticity (str): the shared denticity of all the ligands in the batch.
        optimized (bool, optional): whether the TMCs are made with optimized ligands. Defaults to False.
        df_smiles (pd.DataFrame, optional): reference for decoded SMILES. Defaults to None.

    Returns:
        The extended dataframe.
    """
    if denticity == "monodentate":
        num_ligands = 4
    elif denticity == "bidentate":
        num_ligands = 2

    if optimized:
        ref_smiles_column = "Optimized encoded SMILES"
        ref_xyz_column = "Optimized XYZ"
        ref_smiles_are_encoded = True
    else:
        if denticity == "monodentate":
            ref_smiles_column = "Labeled SMILES"
        elif denticity == "bidentate":
            ref_smiles_column = "Decoded SMILES"
        ref_xyz_column = "XYZ"
        ref_smiles_are_encoded = False

    decoded_smilesx = []
    decoded_mols = []
    log_ps = []
    g_params = []
    for i, row in df.iterrows():
        if df_smiles is not None:
            ligand_id = row["Ligand ID"]
            df_temp = df_smiles[df_smiles["Ligand ID"] == ligand_id]
            assert df_temp.shape[0] == 1
            decoded_smiles = df_temp["Decoded SMILES"].values[0]
            decoded_mol = Chem.MolFromSmiles(decoded_smiles)
        else:
            smiles = row[ref_smiles_column]
            if ref_smiles_are_encoded:
                encoded_mol = Chem.MolFromSmiles(smiles)
                if denticity == "monodentate":
                    decoded_mol, _ = process_substitute_attachment_points(encoded_mol)
                elif denticity == "bidentate":
                    decoded_mol, _ = process_substitute_attachment_points_bidentate(
                        encoded_mol
                    )
                decoded_smiles = Chem.MolToSmiles(decoded_mol)
            else:
                decoded_smiles = smiles
                decoded_mol = Chem.MolFromSmiles(decoded_smiles)
        decoded_smilesx.append(decoded_smiles)
        decoded_mols.append(decoded_mol)
        if decoded_mol is None:
            print(f"{i}: could not make a Mol object from {decoded_smiles}.")
            log_p = None
        else:
            log_p = rdMolDescriptors.CalcCrippenDescriptors(decoded_mol)[0]
        log_ps.append(log_p)

        xyz = row[ref_xyz_column]
        g_param = get_g_parameter(xyz, num_ligands)
        g_params.append(g_param)

    if optimized:
        df["Optimized log P"] = log_ps
        df["Optimized G parameter"] = g_params
    else:
        df["log P"] = log_ps
        df["G parameter"] = g_params

    return df


def batch_parse_orca_logs(
    input_dir: Path,
    denticity: str,
    columns: list,
    optimized: bool = False,
    df_smiles: pd.DataFrame = None,
) -> pd.DataFrame:
    """Parse from a directory a batch of ORCA logfiles named according to the
    respective homoleptic TMCs.

    Args:
        input_dir (Path): directory to ORCA logs of normally terminating TMC calculations.
        denticity (str): the shared denticity of all the ligands in the batch.
        columns (list): the column names to include in the output dataframe.
        optimized (bool): whether the TMCs are made with optimized ligands. Defaults to False.
        df_smiles (pd.DataFrame): reference for decoded SMILES.

    Returns:
        The output dataframe with the parsed TMC labels and optimized geometry.
    """
    files = os.listdir(input_dir)
    files = [f for f in files if os.path.isfile(os.path.join(input_dir, f))]
    if denticity == "monodentate":
        # E.g.: uncond_mono-min15k-356.out
        filename_pattern = "{ligand_id}.out"
    elif denticity == "bidentate":
        # E.g.: uncond_bi-min10k-9850-trans.out
        filename_pattern = (
            "{ligand_series_part1}-{ligand_series_part2}-{n}-{isomer}.out"
        )

    data_dict = {key: [] for key in columns}

    for i, f in enumerate(files):
        if i % 1000 == 0:
            print(f"Starting {i+1} of {len(files)} files.")
        parsed_filename = parse(filename_pattern, f)
        if denticity == "monodentate":
            ligand_id = parsed_filename["ligand_id"]
            isomer = None
        elif denticity == "bidentate":
            ligand_series_part1 = parsed_filename["ligand_series_part1"]
            ligand_series_part2 = parsed_filename["ligand_series_part2"]
            isomer = parsed_filename["isomer"]
            n = parsed_filename["n"]
            ligand_id = f"{ligand_series_part1}-{ligand_series_part2}-{n}"
        filepath = os.path.join(input_dir, f)
        with open(filepath, "r") as fh:
            lines = fh.readlines()
        res = get_orca_results(
            lines,
            properties=[
                "electronic_energy",
                "opt_structure",
                "homo_lumo_energy",
                "cm5",
            ],
        )
        assert res["cm5"]["ATOM"][0] == "Ir"
        metal_charge = res["cm5"]["QCM5"][0]
        homo_lumo_gap = res["homo_lumo_energy"]
        atoms, coords = res["opt_structure"]
        electronic_single_point_energy = res["electronic_energy"]

        if optimized:
            # Eh_to_eV
            data_dict["Optimized ligand ID"].append(ligand_id)
            if isomer is not None:
                data_dict["Optimized isomer"].append(isomer)
            data_dict["Optimized HOMO-LUMO gap (eV)"].append(
                float(homo_lumo_gap) * Eh_to_eV
            )
            data_dict["Optimized metal center charge"].append(metal_charge)
            name = str(f)[:-4]
            data_dict["Optimized XYZ"].append(repackage_xyz(name, atoms, coords))
            data_dict["Optimized final single point energy (Eh)"].append(
                electronic_single_point_energy
            )
        else:
            data_dict["Ligand ID"].append(ligand_id)
            if isomer is not None:
                data_dict["Isomer"].append(isomer)
            data_dict["HOMO-LUMO gap (Eh)"].append(float(homo_lumo_gap))
            data_dict["Metal center charge"].append(metal_charge)
            name = str(f)[:-4]
            data_dict["XYZ"].append(repackage_xyz(name, atoms, coords))
            data_dict["Final single point energy (Eh)"].append(
                electronic_single_point_energy
            )
        df = pd.DataFrame.from_dict(data_dict)
    df = df[columns]
    return df


def combine_optimized_dataframes(
    df_output: pd.DataFrame,
    df_optimized: pd.DataFrame,
    parsing_columns: list,
    columns: list,
) -> pd.DataFrame:
    """Combine the parsed labels from ORCA files with the reference data on the
    respective optimized ligands.

    Args:
        df_output (pd.DataFrame): the results from parsing ORCA files.
        df_optimized (pd.DataFrame): the reference file after optimizing ligands.
        parsing_columns (list): the columns used in parsing ORCA files
        columns (list): the complete set of columns.

    Returns:
        The combined dataframe.
    """
    data_dict = {key: [] for key in columns}
    for i, row in df_output.iterrows():
        optimized_ligand_id = row["Optimized ligand ID"]
        df_temp = df_optimized[
            df_optimized["Optimized ligand ID"] == optimized_ligand_id
        ]
        if "Optimized isomer" in row:
            optimized_isomer = row["Optimized isomer"]
            df_temp = df_temp[df_temp["Optimized isomer"] == optimized_isomer]
        else:
            optimized_isomer = None
        for col in columns:
            if col in parsing_columns:
                data_dict[col].append(row[col])
            else:
                data_dict[col].append(df_temp[col].values[0])

    return pd.DataFrame(data=data_dict, columns=columns)


if __name__ == "__main__":
    args = parse_args()
    denticity = args.denticity
    xtent = args.xtent
    input_dir = os.path.abspath(args.input_dir)
    if args.type == "unoptimized":
        optimized = False
    elif args.type == "optimized":
        optimized = True
        # Required iff optimized:
        reference_optimization_path = os.path.abspath(args.reference_optimization)
        df_optimized = pd.read_csv(reference_optimization_path)
    output_csv_path = os.path.abspath(args.output_file)
    output_xyz_path = os.path.abspath(str(args.output_file)[:-4] + ".xyz")
    df_smiles = None

    if not optimized:
        if denticity == "monodentate":
            columns = [
                "Ligand ID",
                "HOMO-LUMO gap (Eh)",
                "Metal center charge",
                "XYZ",
                "Final single point energy (Eh)",
                "log P",
                "G parameter",
            ]
            parsing_columns = [
                "Ligand ID",
                "HOMO-LUMO gap (Eh)",
                "Metal center charge",
                "XYZ",
                "Final single point energy (Eh)",
            ]
        elif denticity == "bidentate":
            columns = [
                "Ligand ID",
                "Isomer",
                "HOMO-LUMO gap (Eh)",
                "Metal center charge",
                "XYZ",
                "Final single point energy (Eh)",
                "log P",
                "G parameter",
            ]
            parsing_columns = [
                "Ligand ID",
                "Isomer",
                "HOMO-LUMO gap (Eh)",
                "Metal center charge",
                "XYZ",
                "Final single point energy (Eh)",
            ]
        reference_smiles_path = os.path.abspath(args.smiles_reference)
        df_smiles = pd.read_csv(reference_smiles_path)

    elif optimized:
        if denticity == "monodentate":
            columns = [
                "Label ID",
                "Original ligand ID",
                "Optimized ligand ID",
                "Original decoded SMILES",
                "Sampling region",
                "Optimization objective",
                "Minimization",
                "Original encoded SMILES",
                "Optimized encoded SMILES",
                "Tanimoto similarity",
                "Original metal center charge",
                "Optimized metal center charge",
                "Original HOMO-LUMO gap (Eh)",
                "Original HOMO-LUMO gap (eV)",
                "Optimized HOMO-LUMO gap (eV)",
                "Original XYZ",
                "Optimized XYZ",
                "Original final single point energy (Eh)",
                "Optimized final single point energy (Eh)",
                "Original log P",
                "Optimized log P",
                "Original G parameter",
                "Optimized G parameter",
            ]
            parsing_columns = [
                "Optimized ligand ID",
                "Optimized metal center charge",
                "Optimized HOMO-LUMO gap (eV)",
                "Optimized XYZ",
                "Optimized final single point energy (Eh)",
            ]
        elif denticity == "bidentate":
            columns = [
                "Label ID",
                "Original ligand ID",
                "Original isomer",
                "Optimized ligand ID",
                "Optimized isomer",
                "Original decoded SMILES",
                "Sampling region",
                "Optimization objective",
                "Minimization",
                "Original encoded SMILES",
                "Optimized encoded SMILES",
                "Tanimoto similarity",
                "Original metal center charge",
                "Optimized metal center charge",
                "Original HOMO-LUMO gap (Eh)",
                "Original HOMO-LUMO gap (eV)",
                "Optimized HOMO-LUMO gap (eV)",
                "Original XYZ",
                "Optimized XYZ",
                "Original final single point energy (Eh)",
                "Optimized final single point energy (Eh)",
                "Original log P",
                "Optimized log P",
                "Original G parameter",
                "Optimized G parameter",
            ]
            parsing_columns = [
                "Optimized ligand ID",
                "Optimized isomer",
                "Optimized metal center charge",
                "Optimized HOMO-LUMO gap (eV)",
                "Optimized XYZ",
                "Optimized final single point energy (Eh)",
            ]

    df_output = batch_parse_orca_logs(
        input_dir, denticity, parsing_columns, optimized, df_smiles
    )
    df_output = add_ligand_descriptors(df_output, denticity, optimized, df_smiles)
    if optimized:
        xyzs = df_output["Optimized XYZ"].values.tolist()
        df_output = combine_optimized_dataframes(
            df_output, df_optimized, parsing_columns, columns
        )
    elif not optimized:
        xyzs = df_output["XYZ"].values.tolist()

    if xtent == "full":
        df_output.to_csv(output_csv_path, index=None)
        with open(output_xyz_path, "w") as f_out:
            for xyz in xyzs:
                lines = xyz.split("\n")
                for line in lines:
                    f_out.write(line + "\n")
                f_out.write("\n")
    elif xtent == "test":
        df_expect = pd.read_csv(output_csv_path)
        if optimized:
            if denticity == "monodentate":
                df_output = df_output.sort_values(
                    "Optimized ligand ID", axis=0, ascending=True, ignore_index=True
                )
                df_expect = df_expect.sort_values(
                    "Optimized ligand ID", axis=0, ascending=True, ignore_index=True
                )
            elif denticity == "bidentate":
                df_output = df_output.sort_values(
                    ["Optimized ligand ID", "Optimized isomer"],
                    axis=0,
                    ascending=[True, True],
                    ignore_index=True,
                )
                df_expect = df_expect.sort_values(
                    ["Optimized ligand ID", "Optimized isomer"],
                    axis=0,
                    ascending=[True, True],
                    ignore_index=True,
                )

        row_accuracy = compare_dataframes(
            df_output,
            df_expect,
            columns=columns,
        )
        print("")
        if row_accuracy == 1.0:
            print(
                f"Test run perfectly reproduced the pre-existing labeled TMCs file: {output_csv_path}."
            )
        else:
            print(
                f"Test run failed to perfectly reproduce the pre-existing labeled TMCs file: {output_csv_path}."
            )
            print(
                f"Reproduction accuracy of labeled TMCs by matching rows over total rows: {row_accuracy}.\n"
            )

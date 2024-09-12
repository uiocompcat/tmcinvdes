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
                                    -x test
"""

import argparse
import os
from pathlib import Path

import pandas as pd
from parse import parse

from tmcinvdes.quantum_chemistry.orca_parsers import get_orca_results, repackage_xyz
from tmcinvdes.utils import compare_dataframes


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
        default="../../tmp_results/tmp_orca_out/",
        required=True,
        help="Input directory with only successful ORCA .out files, each named after the TMC ID",
    )
    parser.add_argument(
        "--output_file",
        "-o",
        type=Path,
        default="../../tmp_results/uncond_bi-min10k-labeled.csv",
        required=True,
        help="Output filename for .csv, same prefix used for .xyz files.",
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


def batch_parse_orca_logs(
    input_dir: Path, denticity: str, columns: list
) -> pd.DataFrame:
    """Parse from a directory a batch of ORCA logfiles named according to the
    respective TMCs.

    Args:
        input_dir (Path): directory to ORCA logs of normally terminating TMC calculations.
        denticity (str): the shared denticity of all the TMCs in the batch.
        columns (list): the column names to include in the output dataframe.

    Returns:
        The output dataframe with the parsed TMC labels and optimized geometry.
    """
    files = os.listdir(input_dir)
    files = [f for f in files if os.path.isfile(os.path.join(input_dir, f))]
    if denticity == "monodentate":
        filename_pattern = "{ligand_id}.out"  # E.g.: uncond_mono-min15k-356.out
    elif denticity == "bidentate":
        filename_pattern = "{ligand_series_part1}-{ligand_series_part2}-{n}-{isomer}.out"  # E.g.: uncond_bi-min10k-9850-trans.out
        primary_columns = [
            "Ligand ID",
            "Isomer",
            "HOMO-LUMO gap (Eh)",
            "Metal center charge",
            "XYZ",
        ]

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
        data_dict["Ligand ID"].append(ligand_id)
        if isomer is not None:
            data_dict["Isomer"].append(isomer)
        data_dict["HOMO-LUMO gap (Eh)"].append(float(homo_lumo_gap))
        data_dict["Metal center charge"].append(metal_charge)
        atoms, coords = res["opt_structure"]
        name = str(f)[:-4]
        data_dict["XYZ"].append(repackage_xyz(name, atoms, coords))
    df = pd.DataFrame.from_dict(data_dict)
    df = df[primary_columns]
    return df


if __name__ == "__main__":
    args = parse_args()
    denticity = args.denticity
    xtent = args.xtent
    input_dir = os.path.abspath(args.input_dir)
    output_csv_path = os.path.abspath(args.output_file)
    output_xyz_path = os.path.abspath(str(args.output_file)[:-4] + ".xyz")

    if denticity == "monodentate":
        columns = ["Ligand ID", "HOMO-LUMO gap (Eh)", "Metal center charge", "XYZ"]
    elif denticity == "bidentate":
        columns = [
            "Ligand ID",
            "Isomer",
            "HOMO-LUMO gap (Eh)",
            "Metal center charge",
            "XYZ",
        ]

    df_output = batch_parse_orca_logs(input_dir, denticity, columns)
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

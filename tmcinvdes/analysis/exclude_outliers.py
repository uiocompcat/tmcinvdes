"""Apply IsolationForest in order to separate labeled TMC data into outliers to
be excluded and the remaining data points to be included in subsequent training
of conditional generative models.

Usage: (from root directory of this repository)

python -m tmcinvdes.analysis.exclude_outliers -d bidentate \
                                              -i datasets/07_uncond-labeled \
                                              -o datasets/08a_uncond-excluded \
                                              -r datasets/08b_uncond-included \
                                              -x test
"""

import argparse
import os
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd
from scipy import constants
from sklearn.ensemble import IsolationForest

from tmcinvdes.utils import compare_dataframes

eV_per_Eh = 1 / constants.physical_constants["electron volt-hartree relationship"][0]
Eh_per_eV = constants.physical_constants["electron volt-hartree relationship"][0]


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
        help="Input directory with labeled TMCs.",
    )
    parser.add_argument(
        "--output_dir",
        "-o",
        type=Path,
        required=True,
        help="Output directory for the outlier data to be excluded.",
    )
    parser.add_argument(
        "--remainder_dir",
        "-r",
        type=Path,
        required=True,
        help="Output directory for the remainder data to be included.",
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


def remove_outliers(df: pd.DataFrame, contamination: float) -> Tuple[pd.DataFrame]:
    """Remove outliers based on a contamination threshold.

    Args:
        df (pd.DataFrame): the dataframe containing DFT-optimized and -labeled TMC data instances.
        contamination (float): the proportion of data instances to exclude, \in [0.0, 1.0].

    Returns:
        A tuple containing two dataframes:
        First, the remaining (1.0-contamination) data, to be included in subsequent model training,
        after the most extreme instances are excluded.
        Second, the data instances considered as outliers, to be excluded.
    """
    # df["HOMO-LUMO gap (eV)"] = df["HOMO-LUMO gap (Eh)"]*eV_per_Eh
    data_points = df[["Metal center charge", "HOMO-LUMO gap (Eh)"]].to_numpy()
    df_indices = df.index.values
    df_indices = df_indices.reshape(-1, 1)

    # Contamination is percentage of outliers to remove
    iforest = IsolationForest(n_estimators=300, contamination=contamination)
    iforest.fit(data_points)
    iforest_preds = iforest.predict(data_points)
    outlier_indices = np.argwhere(iforest_preds == -1)
    outlier_indices = outlier_indices.flatten()

    outlier_indices = np.append(
        np.zeros(0), df_indices[outlier_indices].flatten(), axis=0
    )

    df_excl = df.iloc[outlier_indices]
    df_incl = df.drop(outlier_indices)
    return df_incl, df_excl


if __name__ == "__main__":
    args = parse_args()
    denticity = args.denticity
    xtent = args.xtent
    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output_dir)
    remainder_dir = os.path.abspath(args.remainder_dir)

    # Modify parameters based on denticity:
    if denticity == "monodentate":
        input_file = os.path.join(input_dir, "uncond_mono-min15k-labeled.csv")
        output_file = os.path.join(
            output_dir, "uncond_mono-min15k-labeled-excluded.csv"
        )
        remainder_file = os.path.join(
            remainder_dir, "uncond_mono-min15k-labeled-included.csv"
        )
        contamination = 0.01
        columns = [
            "Ligand ID",
            "Label ID",
            "Labeled SMILES",
            "HOMO-LUMO gap (Eh)",
            "Metal center charge",
            "XYZ",
        ]
    elif denticity == "bidentate":
        input_file = os.path.join(input_dir, "uncond_bi-min10k-labeled.csv")
        output_file = os.path.join(output_dir, "uncond_bi-min10k-labeled-excluded.csv")
        remainder_file = os.path.join(
            remainder_dir, "uncond_bi-min10k-labeled-included.csv"
        )
        contamination = 0.015
        columns = [
            "Ligand ID",
            "Isomer",
            "Label ID",
            "HOMO-LUMO gap (Eh)",
            "Metal center charge",
            "XYZ",
        ]

    df_input = pd.read_csv(input_file)

    # Modify process based on xtent:
    if xtent == "test":
        df_output_incl, df_output_excl = remove_outliers(
            df_input, contamination=contamination
        )
        df_expect_excl = pd.read_csv(output_file)
        df_expect_incl = pd.read_csv(remainder_file)

        row_accuracy_excl = compare_dataframes(
            df_output_excl, df_expect_excl, columns=columns
        )
        if row_accuracy_excl == 1.0:
            print(
                f"Test run perfectly reproduced the pre-existing excluded outliers file: {output_file}."
            )
        else:
            print(
                f"Test run failed to perfectly reproduce the pre-existing excluded outliers file: {output_file}."
            )
            print(
                f"Reproduction accuracy of excluded outliers by matching rows over total rows: {row_accuracy_excl}."
            )
        print("")
        row_accuracy_incl = compare_dataframes(
            df_output_incl, df_expect_incl, columns=columns
        )
        if row_accuracy_incl == 1.0:
            print(
                f"Test run perfectly reproduced the pre-existing included remainders file: {remainder_file}."
            )
        else:
            print(
                f"Test run failed to perfectly reproduce the pre-existing included remainders file: {remainder_file}."
            )
            print(
                f"Reproduction accuracy of included remainders by matching rows over total rows: {row_accuracy_incl}."
            )
    elif xtent == "full":
        df_output_incl, df_output_excl = remove_outliers(
            df_input, contamination=contamination
        )
        df_output_excl.to_csv(output_file, index=False, header=True)
        df_output_incl.to_csv(remainder_file, index=False, header=True)

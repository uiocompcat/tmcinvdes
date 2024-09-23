"""Sample labeled homoleptic TMCs to select their ligands for subsequent
optimization with conditional models.

This script also produces plots with sampled TMCs/ligands in the foreground and the mass of labeled
instances in the background.

Specifically, this script assumes as input labeled instances that remain included after outliers
were excluded in the previous stage.

The symbols $\epsilon$ and $q_\text{Ir}$ are used in the plots and discussions of the associated
2D plane to represent, respectively, the HOMO-LUMO gap and the metal center charge of the TMCs that
were labeled using DFT calculations.

Usage: (from root directory of repository)

python -m tmcinvdes.ligand_generation.sample_ligands_to_optimize -d bidentate \
                                                            -i datasets/08b_uncond-included/ \
                                                            -o datasets/10_cond-sampled_from_8b/ \
                                                            -x test
"""

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import constants

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
        "--configuration",
        "-c",
        choices=["cis", "trans", "cistrans", ""],
        default="",
        help="Select one type of isomerism supported, or leave as empty string.",
    )
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
        help="Directory with .CSV data files of labeled, included TMCs.",
    )
    parser.add_argument(
        "--output_dir",
        "-o",
        type=Path,
        required=True,
        help="Directory for .CSV data files of sampled ligands to optimize, and for plot images.",
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


def sample_monodentate(df_input: pd.DataFrame) -> tuple[plt.figure, pd.DataFrame]:
    """Sample monodentate ligands from labeled, included TMCs based on HOMO-
    LUMO gap (eV) and central metal charge.

    Args:
        df_input (pd.DataFrame): the complete set of labeled, included homoleptic TMCs with
        monodentate ligands.

    Returns:
        The DataFrame of the selected instances from the input dataframe.
        The Figure object which may or may not be written to file based on xtent.
    """
    range_homo_lumo_gap = [1, 6]
    range_metal_charge = [-0.2, 1]

    min_homo_lumo_gap = [1, 2]
    max_homo_lumo_gap = [4.5, 6]

    min_homo_lumo_gap_alternate = [1, 3]

    min_metal_charge = [-0.2, 0]
    max_metal_charge = [0.45, 0.7]

    center_homo_lumo_gap_upper = [3.51, 4]
    center_homo_lumo_gap_lower = [3, 3.49]

    center_metal_charge_upper = [0.201, 0.25]
    center_metal_charge_lower = [0.15, 0.199]

    mapper = {
        "max_metal_charge": (max_metal_charge, range_homo_lumo_gap),
        "min_metal_charge": (min_metal_charge, range_homo_lumo_gap),
        "max_homo_lumo_gap": (range_metal_charge, max_homo_lumo_gap),
        "min_homo_lumo_gap": (range_metal_charge, min_homo_lumo_gap),
        "center_hi_homo_lumo_gap_hi_metal_charge": (
            center_metal_charge_upper,
            center_homo_lumo_gap_upper,
        ),
        "center_lo_homo_lumo_gap_hi_metal_charge": (
            center_metal_charge_upper,
            center_homo_lumo_gap_lower,
        ),
        "center_lo_homo_lumo_gap_lo_metal_charge": (
            center_metal_charge_lower,
            center_homo_lumo_gap_lower,
        ),
        "center_hi_homo_lumo_gap_lo_metal_charge": (
            center_metal_charge_lower,
            center_homo_lumo_gap_upper,
        ),
    }
    alternate_mapper = {
        "min_metal_charge_min_homo_lumo_gap_alternate": (
            min_metal_charge,
            min_homo_lumo_gap_alternate,
        ),
        "max_metal_charge_min_homo_lumo_gap_alternate": (
            max_metal_charge,
            min_homo_lumo_gap_alternate,
        ),
        "center_alternate": (),
    }

    fig, ax = plt.subplots(figsize=(12, 6))

    # Prepare labeled, included instances for 2D plotting.
    metal_charges = df_input["Metal center charge"].to_numpy()
    df_input["HOMO-LUMO gap (eV)"] = df_input["HOMO-LUMO gap (Eh)"] * eV_per_Eh
    homo_lumo_gaps = df_input["HOMO-LUMO gap (eV)"].to_numpy()
    ax.scatter(metal_charges, homo_lumo_gaps, label="Raw data", color="k")

    df_total = pd.DataFrame([], columns=df_input.columns)

    for i, (k, v) in enumerate(mapper.items()):
        mod = 0
        limits = v
        print(k, limits)
        if k == "max_homo_lumo_gap":
            mod = 0.3
        snip = df_input[
            df_input["Metal center charge"].between(
                limits[0][0], limits[0][1], inclusive="both"
            )
        ][
            df_input["HOMO-LUMO gap (eV)"].between(
                limits[1][0] + mod, limits[1][1], inclusive="both"
            )
        ]

        sample_size = 20
        try:
            sample = snip.sample(sample_size, random_state=1)
        except Exception:
            print(f"Warning for {k}")
            limits = alternate_mapper[f"{k}_alternate"]
            snip = df_input[
                df_input["Metal center charge"].between(
                    limits[0][0], limits[0][1], inclusive="both"
                )
            ][
                df_input["HOMO-LUMO gap (eV)"].between(
                    limits[1][0], limits[1][1], inclusive="both"
                )
            ]
            sample = snip.sample(sample_size, random_state=1)

        if k.startswith("center"):
            sampling_region = "Center"
        elif k == "max_homo_lumo_gap":
            sampling_region = "High HOMO-LUMO gap (eV)"
        elif k == "min_homo_lumo_gap":
            sampling_region = "Low HOMO-LUMO gap (eV)"
        elif k == "max_metal_charge":
            sampling_region = "High metal center charge"
        elif k == "min_metal_charge":
            sampling_region = "Low metal center charge"

        sample["Sampling region"] = [f"{sampling_region}"] * len(sample)
        df_total = pd.concat([df_total, sample])

    for sampling_region in [
        "High metal center charge",
        "Low metal center charge",
        "High HOMO-LUMO gap (eV)",
        "Low HOMO-LUMO gap (eV)",
        "Center",
    ]:
        df_sample = df_total[df_total["Sampling region"] == sampling_region]
        ax.scatter(
            df_sample["Metal center charge"],
            df_sample["HOMO-LUMO gap (eV)"],
            label=f"{sampling_region}",
            s=80,
        )

    ax.set_xlim([-0.5, 1])
    ax.set(
        xlabel="$q_\\text{Ir}$: Metal center charge",
        ylabel="$\epsilon$: HOMO-LUMO gap (eV)",
        title="",
    )
    ax.legend()
    plt.show()

    return df_total, fig


def set_boundaries_quantiles(
    df: pd.DataFrame, quantile_points: list
) -> tuple[dict, list]:
    """Set the boundary boxes of the regions in q-epsilon space to sample
    ligands from.

    Args:
        df (pd.DataFrame): dataframe with HOMO-LUMO gap and metal center charge data.
        quantile_points (list): relative points to calculate boundary boxes based on data.

    Returns:
        The sampling region boundaries as a dict.
        The extremal values of the data points as a list of lists of floats.
    """
    homo_lumo_gaps = df["HOMO-LUMO gap (eV)"]
    metal_charges = df["Metal center charge"]
    extremes = [
        [metal_charges.min(), metal_charges.max()],
        [homo_lumo_gaps.min(), homo_lumo_gaps.max()],
    ]

    homo_lumo_gap_quantiles = np.quantile(homo_lumo_gaps, quantile_points)
    min_homo_lumo_gap = homo_lumo_gaps.min()
    max_homo_lumo_gap = homo_lumo_gaps.max()
    print(
        f"Smallest to greatest HOMO-LUMO gap: [{min_homo_lumo_gap}, {max_homo_lumo_gap}]"
    )
    min_homo_lumo_gap = homo_lumo_gap_quantiles[0]
    max_homo_lumo_gap = homo_lumo_gap_quantiles[4]
    quant_1_homo_lumo_gap = homo_lumo_gap_quantiles[1]
    quant_3_homo_lumo_gap = homo_lumo_gap_quantiles[3]
    mid_homo_lumo_gap = homo_lumo_gap_quantiles[2]
    inner_radius_homo_lumo_gap = (
        max(
            [
                mid_homo_lumo_gap - quant_1_homo_lumo_gap,
                quant_3_homo_lumo_gap - mid_homo_lumo_gap,
            ]
        )
        / 4
    )

    metal_charge_quantiles = np.quantile(metal_charges, quantile_points)
    min_metal_charge = metal_charges.min()
    max_metal_charge = metal_charges.max()
    print(
        f"Smallest to greatest metal center charge: [{min_metal_charge}, {max_metal_charge}]"
    )
    min_metal_charge = metal_charge_quantiles[0]
    max_metal_charge = metal_charge_quantiles[4]
    quant_1_metal_charge = metal_charge_quantiles[1]
    quant_3_metal_charge = metal_charge_quantiles[3]
    mid_metal_charge = metal_charge_quantiles[2]
    inner_radius_metal_charge = (
        max(
            [
                mid_metal_charge - quant_1_metal_charge,
                quant_3_metal_charge - mid_metal_charge,
            ]
        )
        / 4
    )

    ranges_low_homo_lumo_gaps = [min_homo_lumo_gap, quant_1_homo_lumo_gap]
    ranges_high_homo_lumo_gaps = [quant_3_homo_lumo_gap, max_homo_lumo_gap]
    ranges_inner_homo_lumo_gaps_lower = [
        mid_homo_lumo_gap - inner_radius_homo_lumo_gap,
        mid_homo_lumo_gap,
    ]
    ranges_inner_homo_lumo_gaps_higher = [
        mid_homo_lumo_gap,
        mid_homo_lumo_gap + inner_radius_homo_lumo_gap,
    ]
    ranges_inner_homo_lumo_gaps = [
        quant_1_homo_lumo_gap,
        quant_3_homo_lumo_gap,
    ]

    ranges_low_metal_charges = [min_metal_charge, quant_1_metal_charge]
    ranges_high_metal_charges = [quant_3_metal_charge, max_metal_charge]
    ranges_inner_metal_charges_lower = [
        mid_metal_charge - inner_radius_metal_charge,
        mid_metal_charge,
    ]
    ranges_inner_metal_charges_higher = [
        mid_metal_charge,
        mid_metal_charge + inner_radius_metal_charge,
    ]
    ranges_inner_metal_charges = [
        quant_1_metal_charge,
        quant_3_metal_charge,
    ]

    mapper = {
        "high_metal_charge": (ranges_high_metal_charges, ranges_inner_homo_lumo_gaps),
        "low_metal_charge": (ranges_low_metal_charges, ranges_inner_homo_lumo_gaps),
        "high_homo_lumo_gap": (ranges_inner_metal_charges, ranges_high_homo_lumo_gaps),
        "low_homo_lumo_gap": (ranges_inner_metal_charges, ranges_low_homo_lumo_gaps),
        "center_hi_homo_lumo_gap_hi_metal_charge": (
            ranges_inner_metal_charges_higher,
            ranges_inner_homo_lumo_gaps_higher,
        ),
        "center_lo_homo_lumo_gap_hi_metal_charge": (
            ranges_inner_metal_charges_higher,
            ranges_inner_homo_lumo_gaps_lower,
        ),
        "center_lo_homo_lumo_gap_lo_metal_charge": (
            ranges_inner_metal_charges_lower,
            ranges_inner_homo_lumo_gaps_lower,
        ),
        "center_hi_homo_lumo_gap_lo_metal_charge": (
            ranges_inner_metal_charges_lower,
            ranges_inner_homo_lumo_gaps_higher,
        ),
    }

    return mapper, extremes


def sample_bidentate(df_input: pd.DataFrame) -> tuple[plt.figure, pd.DataFrame]:
    """Sample bidentate ligands from labeled, included TMCs based on HOMO-LUMO
    gap (eV) and central metal charge.

    Args:
        df_input (pd.DataFrame): the complete set of labeled, included homoleptic TMCs with
        bidentate ligands.

    Returns:
        The DataFrame of the selected instances from the input dataframe.
        The Figure object which may or may not be written to file based on xtent.
    """
    # Number of samples to take in each region:
    n_samples = 20
    # Small numbers to set the bounding boxes of sampling regions as quantiles of the dataset:
    m = 0.0
    n = 0.1  # m=0.0 and n=0.1 represents the lower and upper 10% of a given range.
    quantile_points = [
        m,
        m + n,
        0.5,
        1.0 - (n + m),
        1.0 - m,
    ]

    df_input["HOMO-LUMO gap (eV)"] = df_input["HOMO-LUMO gap (Eh)"] * eV_per_Eh
    mapper, extremes = set_boundaries_quantiles(df_input, quantile_points)
    homo_lumo_gaps = df_input["HOMO-LUMO gap (eV)"].to_numpy()
    metal_charges = df_input["Metal center charge"].to_numpy()
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.scatter(metal_charges, homo_lumo_gaps, label="Raw data", color="k")
    df_total = pd.DataFrame([], columns=df_input.columns)

    for key, limits in mapper.items():
        df_region = df_input[
            df_input["Metal center charge"].between(
                limits[0][0], limits[0][1], inclusive="both"
            )
            & df_input["HOMO-LUMO gap (eV)"].between(
                limits[1][0], limits[1][1], inclusive="both"
            )
        ]

        df_sample = df_region.sample(n_samples, random_state=1)
        df_sample = df_sample.reset_index(drop=True)

        if key.startswith("center"):
            sampling_region = "Center"
        elif key == "high_homo_lumo_gap":
            sampling_region = "High HOMO-LUMO gap (eV)"
        elif key == "low_homo_lumo_gap":
            sampling_region = "Low HOMO-LUMO gap (eV)"
        elif key == "high_metal_charge":
            sampling_region = "High metal center charge"
        elif key == "low_metal_charge":
            sampling_region = "Low metal center charge"
        df_sample["Sampling region"] = [f"{sampling_region}"] * len(df_sample)
        df_total = pd.concat([df_total, df_sample])
    for sampling_region in df_total["Sampling region"].unique():
        df_sample = df_total[df_total["Sampling region"] == sampling_region]
        ax.scatter(
            df_sample["Metal center charge"],
            df_sample["HOMO-LUMO gap (eV)"],
            label=f"{sampling_region}",
            s=80,
        )
    min_metal_charge = extremes[0][0]
    max_metal_charge = extremes[0][1]
    delta_metal_charge = (max_metal_charge - min_metal_charge) * 0.05
    ax.set_xlim(
        [min_metal_charge - delta_metal_charge, max_metal_charge + delta_metal_charge]
    )
    ax.set(
        xlabel="$q_\\text{Ir}$: Metal center charge",
        ylabel="$\epsilon$: HOMO-LUMO gap (eV)",
        title="",
    )

    ax.legend()
    plt.show()

    df_total = df_total.drop_duplicates(subset=["Ligand ID", "Isomer"])
    return df_total, fig


if __name__ == "__main__":
    args = parse_args()
    isomer_config = args.configuration
    denticity = args.denticity
    xtent = args.xtent
    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output_dir)

    # Modify parameters based on denticity:
    if denticity == "monodentate":
        input_file = os.path.join(input_dir, "uncond_mono-min15k-labeled-included.csv")
        output_csv_file = os.path.join(
            output_dir, "uncond_mono-min15k-labeled-included-sampled_for_cond_mono.csv"
        )
        output_png_file = os.path.join(
            output_dir, "uncond_mono-min15k-labeled-included-sampled_for_cond_mono.png"
        )
        columns = [
            "Ligand ID",
            "Label ID",
            "Labeled SMILES",
            "HOMO-LUMO gap (Eh)",
            "Metal center charge",
            "XYZ",
        ]
    elif denticity == "bidentate":
        if isomer_config == "":
            isomer_config = "cis_trans"
        elif isomer_config == "cistrans":
            isomer_config = "cis_trans"
        input_file = os.path.join(input_dir, "uncond_bi-min10k-labeled-included.csv")
        output_csv_file = os.path.join(
            output_dir,
            f"uncond_bi-min10k-labeled-included-sampled_for_cond_bi_{isomer_config}.csv",
        )
        output_png_file = os.path.join(
            output_dir,
            f"uncond_bi-min10k-labeled-included-sampled_for_cond_bi_{isomer_config}.png",
        )
        columns = [
            "Ligand ID",
            "Isomer",
            "Label ID",
            "HOMO-LUMO gap (Eh)",
            "Metal center charge",
            "XYZ",
        ]

    df_input = pd.read_csv(input_file)

    # Perform sampling and plotting:
    if denticity == "monodentate":
        df_output, fig = sample_monodentate(df_input)
    elif denticity == "bidentate":
        df_output, fig = sample_bidentate(df_input)

    # Modify process based on xtent:
    if xtent == "test":
        df_expect = pd.read_csv(output_csv_file)

        row_accuracy = compare_dataframes(df_output, df_expect, columns=columns)
        if row_accuracy == 1.0:
            print(
                f"Test run perfectly reproduced the pre-existing sampled ligands file: {output_csv_file}."
            )
        else:
            print(
                f"Test run failed to perfectly reproduce the pre-existing sampled ligands file: {output_csv_file}."
            )
            print(
                f"Reproduction accuracy of sampled ligands by matching rows over total rows: {row_accuracy}."
            )
        print("")
    elif xtent == "full":
        print(df_output.shape)
        df_output.reset_index(drop=True, inplace=True)
        df_output.to_csv(output_csv_file, index=False, header=True)
        fig.savefig(output_png_file, dpi=600)

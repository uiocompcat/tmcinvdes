"""Assembles single-center iridium TMCs from a generated set of ligands.

Can take either monodentate ligands to make [IrL_4]^+ complexes, or bidentate ligands to make
[IrL_2]^+ complexes as both cis and trans isomers.

Usage: (from root directory of repository)

python -m tmcinvdes.structure_generation.assemble_tmcs -d bidentate \
                                                       -i datasets/05_uncond_minXk/uncond_bi-min10k.csv \
                                                       -o datasets/06_uncond-TMC/uncond_bi-min10k-TMC.xyz \
                                                       -x demo
"""

import argparse
import os
import shutil
from datetime import datetime
from pathlib import Path

import pandas as pd

from tmcinvdes.structure_generation.batch_add_ligands import (
    batch_add_ligands,
    batch_add_monodentate_optimized_ligands,
)
from tmcinvdes.structure_generation.molsimplify_tools import run_bash


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
        "--input_file",
        "-i",
        type=Path,
        required=True,
        help="Input file with correct data, describe expected input file here.",
    )
    parser.add_argument(
        "--output_file",
        "-o",
        type=Path,
        required=True,
        help="Output file(s), describe briefly.",
    )
    parser.add_argument(
        "--optimized",
        "-p",
        type=bool,
        default=False,
        help="Indicate whether ligands in input file were optimized.",
    )
    parser.add_argument(
        "--xtent",
        "-x",
        choices=["full", "test", "demo"],
        default="test",
        help="""The extent to which the process should run.
         - In the "full" case, the process runs in its entirety and writes the output results to
           file.
         - In the "test" case, the process runs in its entirety and compares the results to the
           output already existing in the designated output file.
         - In the "demo" case, the process only runs to a limited extent (n ~ 100 instances) and
           if applicable compares the limited results to the corresponding instances already
           existing in the designated output file.""",
    )
    return parser.parse_args(arg_list)


def build_tmc(geometry: str, metal_center: str, ligand_names: list) -> tuple:
    """Build a TMC from ligands in the local ligand library.

    Args:
        geometry (str): The metal complex geometry (li, tpl, sqp, thd, spy, tbp, oct, tpr, pbp, sqap).
        metal_center (str): The metal center element identifier.
        ligand_names (list): List of the ligand names.

    Returns:
        A tuple of the following two variables:
        xyz (str): the XYZ of the built metal complex.
        status (str): the status of the complex generation run.
    """
    # Remove previous (if any) molSimplify operation directory.
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

    output = run_bash(" ".join(["molsimplify", *parameters]))

    # Save the molSimplify structure.
    if output.find("WARNING: Generated complex is not good!") > -1:
        status = "bad"
    else:
        status = "good"
    try:
        with open("Runs/run/run/run.xyz") as fh:
            xyz = fh.read()
    except Exception:
        xyz = None
        status = "bad"
    return xyz, status


def get_xyz_dict(xyzs: list) -> dict:
    """Make a dictionary of XYZs where the comment line is the key and the
    whole XYZ structure (text string) is the value.

    Args:
        xyzs (list): The list of XYZs.

    Returns:
        A dictionary of XYZs.
    """
    xyz_dict = {}
    for xyz in xyzs:
        key = xyz.split("\n")[1]
        xyz_dict[key] = xyz.strip()

    return xyz_dict


def get_elements(xyz: str) -> list:
    """Extract the elements from a given XYZ structure.

    Args:
        xyz (str): The XYZ structure.

    Returns:
        The list of elements.
    """
    return [_.split()[0] for _ in xyz.split("\n")[2:]]


def compare_xyzs(xyzs_a: list, xyzs_b: list) -> bool:
    """Compare two lists of XYZ structures.

    Args:
        xyzs_a (list): The first list of XYZ structures.
        xyzs_b (list): The second list of XYZ structures.

    Returns:
        A boolean flag indicating whether the two lists of XYZ structures are equivalent.
    """
    xyz_dict_a = get_xyz_dict(xyzs_a)
    xyz_dict_b = get_xyz_dict(xyzs_b)

    xyz_dict_a_keys = set(list(xyz_dict_a.keys()))
    xyz_dict_b_keys = set(list(xyz_dict_b.keys()))

    # Check for the same IDs in both lists.
    if xyz_dict_a_keys != xyz_dict_b_keys:
        print("Differing structures detected:")
        print(
            sorted(
                list(
                    xyz_dict_a_keys.difference(xyz_dict_b_keys)
                    | xyz_dict_b_keys.difference(xyz_dict_a_keys)
                )
            )
        )
    keys_accuracy = len(xyz_dict_a_keys.intersection(xyz_dict_b_keys)) / len(
        xyz_dict_a_keys.union(xyz_dict_b_keys)
    )

    # Check for the same elements per XYZ structure.
    elements_matched_keys = set()
    for key in xyz_dict_a_keys.intersection(xyz_dict_b_keys):
        elements_a = get_elements(xyz_dict_a[key])
        elements_b = get_elements(xyz_dict_b[key])

        if elements_a == elements_b:
            elements_matched_keys.add(key)
        else:
            print("Different elements in XYZ detected in structure: " + key)
    elements_accuracy = len(elements_matched_keys) / len(
        xyz_dict_a_keys.intersection(xyz_dict_b_keys)
    )

    return keys_accuracy, elements_accuracy


def update_progress_df_and_file(
    df_progress: pd.DataFrame, df_temp: pd.DataFrame, progress_file: Path
) -> pd.DataFrame:
    """Update both DataFrame in memory and corresponding .csv on disk for each
    TMC assembly attempted.

    Args:
        df_progress (pd.DataFrame): the dataframe representing total progress thus far.
        df_temp (pd.DataFrame): the smaller dataframe of recent progress to append to both file
        and progress dataframe.
        progress_file (Path): the filepath where the progress is being saved.

    Returns:
        pd.DataFrame: the updated in-memory progress dataframe.
    """
    if not os.path.isfile(progress_file):
        df_temp.to_csv(progress_file, index=False, header=True)
        df_progress = df_temp
    else:
        df_temp.to_csv(progress_file, mode="a", index=False, header=False)
        df_progress = pd.concat([df_progress, df_temp])
    return df_progress


if __name__ == "__main__":
    args = parse_args()
    denticity = args.denticity
    xtent = args.xtent
    input_file = os.path.abspath(args.input_file)
    output_file = os.path.abspath(args.output_file)
    df_input = pd.read_csv(input_file)

    print("Adding ligands to local molSimplify library...")
    # Add ligands to molsimplify library.
    if args.optimized:
        tmp_subdir = "12_cond-TMC"
        ligand_id_key = "Optimized ligand ID"
        if denticity == "monodentate":
            batch_add_monodentate_optimized_ligands(df_input, 0)
    else:
        tmp_subdir = "06_uncond-TMC"
        ligand_id_key = "Ligand ID"
        batch_add_ligands(df_input, 0)
    print("Done!")

    # Add file functionality to support tracking progress and resume more easily.
    output_basename = str(os.path.basename(output_file))
    os.makedirs(f"tmp_results/{tmp_subdir}/", exist_ok=True)

    progress_file = (
        f"tmp_results/{tmp_subdir}/" + output_basename[:-4] + f"_{xtent}_progress.csv"
    )
    progress_file = os.path.join(*progress_file.split("/"))
    if os.path.isfile(progress_file) and xtent in ["full", "test"]:
        df_progress = pd.read_csv(progress_file)
    else:
        df_progress = None

    print("Assembling complexes...")
    # Separate parallel branches.
    if denticity == "monodentate":
        if df_progress is None:
            df_progress = pd.DataFrame(
                [], columns=[ligand_id_key, "TMC assembly status", "XYZ"]
            )
            offset = None
        else:
            offset = df_progress.iloc[-1][ligand_id_key]
            print(f"Resuming complex assembly after {offset}.")

        # Modify input based on xtent.
        if xtent == "demo":  # Constrain how many instances to process.
            if args.optimized:
                demo_ligand_ids = [
                    "cond_mono-O6_optimized-uncond_mono-min15k-352",
                    "cond_mono-O4_optimized-uncond_mono-min15k-1464",
                    "cond_mono-O5_optimized-uncond_mono-min15k-4135",
                    "cond_mono-O4_optimized-uncond_mono-min15k-4357",
                    "cond_mono-O3_optimized-uncond_mono-min15k-4842",
                ]
            else:
                demo_ligand_ids = [
                    "uncond_mono-min15k-1",
                    "uncond_mono-min15k-2",
                    "uncond_mono-min15k-3",
                    "uncond_mono-min15k-4",
                    "uncond_mono-min15k-5",
                ]
            df_input = df_input[df_input[ligand_id_key].isin(demo_ligand_ids)]
        if offset:
            xyzs = df_progress["XYZ"].values.tolist()
            tried_ligand_ids = df_progress[ligand_id_key].values.tolist()
            df_input = df_input[~df_input[ligand_id_key].isin(tried_ligand_ids)]

        # Build homoleptic TMCs and record XYZs.
        if offset is None:
            xyzs = []
        for ligand_id in df_input[ligand_id_key].to_numpy():
            xyz, status = build_tmc(
                geometry="sqp",
                metal_center="Ir",
                ligand_names=4 * [ligand_id],
            )
            if xyz is not None:
                xyz_lines = xyz.split("\n")
                xyz_lines[1] = ligand_id
                xyz = "\n".join(xyz_lines)
            df_temp = pd.DataFrame(
                data={
                    ligand_id_key: [ligand_id],
                    "TMC assembly status": [status],
                    "XYZ": [xyz],
                }
            )
            df_progress = update_progress_df_and_file(
                df_progress, df_temp, progress_file
            )

            # Skip to next ligand if build status is bad.
            if status == "bad":
                continue

            # Append XYZs if building succeeded.
            print(f"Assembled homoleptic TMC with ligand {ligand_id}.")
            xyzs.append(xyz)

    elif denticity == "bidentate":
        if df_progress is None:
            df_progress = pd.DataFrame(
                [], columns=[ligand_id_key, "Isomer", "TMC assembly status", "XYZ"]
            )
            offset = None
            offset_isomer = None
        else:
            offset = df_progress.iloc[-1][ligand_id_key]
            offset_isomer = df_progress.iloc[-1]["Isomer"]
            print(f"Resuming complex assembly after {offset}-{offset_isomer}.")
        # Modify input based on xtent.
        if xtent == "demo":  # Constrain how many instances to process.
            if args.optimized:
                demo_ligand_ids = (
                    []
                )  # Bidentate optimized ligands were not part of the preprint v1.0 work.
            else:
                demo_ligand_ids = [
                    "uncond_bi-min10k-2",
                    "uncond_bi-min10k-3",
                    "uncond_bi-min10k-4",
                    "uncond_bi-min10k-6",
                    "uncond_bi-min10k-7",
                ]
            df_input = df_input[df_input[ligand_id_key].isin(demo_ligand_ids)]
        if offset:
            xyzs = df_progress["XYZ"].values.tolist()
            tried_ligand_ids = df_progress[ligand_id_key].values.tolist()
            df_input = df_input[~df_input[ligand_id_key].isin(tried_ligand_ids)]

        # Build homoleptic TMCs (cis and trans) and record XYZ.
        if offset is None:
            xyzs = []
        for ligand_id in df_input[ligand_id_key].to_numpy():
            # Build cis TMC.
            if offset_isomer in [None, "trans"]:
                offset_isomer = "cis"
                cis_xyz, cis_status = build_tmc(
                    geometry="sqp",
                    metal_center="Ir",
                    ligand_names=[ligand_id, ligand_id + "_flipped"],
                )

                if cis_xyz is not None:
                    cis_xyz_lines = cis_xyz.split("\n")
                    cis_xyz_lines[1] = ligand_id + "-cis"
                    cis_xyz = "\n".join(cis_xyz_lines)

            df_temp = pd.DataFrame(
                data={
                    ligand_id_key: [ligand_id],
                    "Isomer": [offset_isomer],
                    "TMC assembly status": [cis_status],
                    "XYZ": [cis_xyz],
                }
            )
            df_progress = update_progress_df_and_file(
                df_progress, df_temp, progress_file
            )

            # Skip to next ligand if build status is bad.
            if cis_status == "bad":
                offset_isomer = "trans"
                df_temp = pd.DataFrame(
                    data={
                        ligand_id_key: [ligand_id],
                        "Isomer": [offset_isomer],
                        "TMC assembly status": ["skipped"],
                        "XYZ": [None],
                    }
                )
                df_progress = update_progress_df_and_file(
                    df_progress, df_temp, progress_file
                )
                continue

            # Build trans TMC.
            if offset_isomer == "cis":
                offset_isomer = "trans"
                trans_xyz, trans_status = build_tmc(
                    geometry="sqp",
                    metal_center="Ir",
                    ligand_names=[ligand_id, ligand_id],
                )

                if trans_xyz is not None:
                    trans_xyz_lines = trans_xyz.split("\n")
                    trans_xyz_lines[1] = ligand_id + "-trans"
                    trans_xyz = "\n".join(trans_xyz_lines)

                df_temp = pd.DataFrame(
                    data={
                        ligand_id_key: [ligand_id],
                        "Isomer": [offset_isomer],
                        "TMC assembly status": [trans_status],
                        "XYZ": [trans_xyz],
                    }
                )
                df_progress = update_progress_df_and_file(
                    df_progress, df_temp, progress_file
                )

                # Skip to next ligand if build status is bad.
                if trans_status == "bad":
                    continue

            # Append XYZs if building succeeded for both cis and trans isomerism.
            print(f"Assembled both cis and trans TMCs with ligand {ligand_id}.")
            xyzs.append(cis_xyz)
            xyzs.append(trans_xyz)

    # Modify output based on xtent.
    if xtent == "full":
        # Write current results to file. This will overwrite current output_file.
        with open(output_file, "w") as fh:
            fh.write("\n\n".join(xyzs))
    elif xtent == "test":
        # Compare the current results with results from disk.
        try:
            with open(output_file, "r") as fh:
                expected_xyzs = fh.read().split("\n\n")
        except Exception:
            pass

        keys_accuracy, elements_accuracy = compare_xyzs(xyzs, expected_xyzs)

        # Since the full process is so extensive, write test results to file.

        with open(
            f"tmp_results/{tmp_subdir}/test-"
            + datetime.now().strftime("%m%d%Y-%H%M%S")
            + ".xyz",
            "w",
        ) as fh:
            fh.write("\n\n".join(xyzs))
        with open(
            f"tmp_results/{tmp_subdir}/test-"
            + datetime.now().strftime("%m%d%Y-%H%M%S")
            + ".txt",
            "w",
        ) as fh:
            fh.write(
                f"Accuracy keys: {keys_accuracy}\nAccuracy elements: {elements_accuracy}"
            )

        if keys_accuracy == 1.0 and elements_accuracy == 1.0:
            print("The process reproduced the expected results entirely!")
        else:
            print(
                f"When testing the reproducibility, the newly assembled TMCs matched (in terms of ligands and isomerism used) the set of previously assembled TMCs with an accuracy of {keys_accuracy:.3f}."
            )
            print(
                f"When testing the reproducibility, of TMCs with matching ligand ID and isomerism present among both newly and previously assembled TMCs, the matched TMCs also matched (in terms of constituent atoms' element, order, and number) with an accuracy of {elements_accuracy:.3f}."
            )
    elif xtent == "demo":
        # Compare the current results with results from disk for only the selected instances.
        try:
            with open(output_file, "r") as fh:
                expected_xyzs = [
                    _
                    for _ in fh.read().split("\n\n")
                    if _.split("\n")[1].replace("-cis", "").replace("-trans", "")
                    in df_input[ligand_id_key].to_numpy().tolist()
                ]
        except Exception:
            pass

        keys_accuracy, elements_accuracy = compare_xyzs(xyzs, expected_xyzs)

        # Write demo reproduction results to file.
        os.makedirs(f"tmp_results/{tmp_subdir}/", exist_ok=True)
        with open(
            f"tmp_results/{tmp_subdir}/demo-"
            + datetime.now().strftime("%m%d%Y-%H%M%S")
            + ".xyz",
            "w",
        ) as fh:
            fh.write("\n\n".join(xyzs))
        with open(
            f"tmp_results/{tmp_subdir}/demo-"
            + datetime.now().strftime("%m%d%Y-%H%M%S")
            + ".txt",
            "w",
        ) as fh:
            fh.write(
                f"Accuracy keys: {keys_accuracy}\nAccuracy elements: {elements_accuracy}"
            )

        if keys_accuracy == 1.0 and elements_accuracy == 1.0:
            print("The process reproduced the expected results entirely!")
        else:
            print(
                f"When testing the reproducibility, the newly assembled TMCs matched (in terms of ligands and isomerism used) the set of previously assembled TMCs with an accuracy of {keys_accuracy:.3f}."
            )
            print(
                f"When testing the reproducibility, of TMCs with matching ligand ID and isomerism present among both newly and previously assembled TMCs, the matched TMCs also matched (in terms of constituent atoms' element, order, and number) with an accuracy of {elements_accuracy:.3f}."
            )

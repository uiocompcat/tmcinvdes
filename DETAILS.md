# Detailed Experimental Workflow

## Overview

Every data workflow stage under [datasets/DETAILS.md](datasets/DETAILS.md) that relies on the code in the present repository can be reproduced by one Python script (`tmcinvdes/.../{process_verb}.py`) being called by two Bash scripts (`bash/{process_verb}-{denticity}-{specifications}.sh`). Reproduction is not always perfect, and further details on this aspect are available [here](REPRODUCIBILITY.md). *(The `-{specifications}` can be used to distinguish Bash files where the Python script is used to produce more than one stage, i.e., for assembling TMC structures, and for parsing ORCA output files to label TMCs/ligands after ORCA calculations.)*

In other words, the Python script structures the general process for both the monodentate and bidentate cases, whereas the Bash scripts represent the nominal run of that process for the monodentate or bidentate case.

## Extent (`--xtent (-x)`)

We use the flag `--xtent` or `-x` to differentiate the extent to which the Python script should run. This generally differentiates between writing results to disk or testing in-memory reproduced results.

- In the `full` case, the process runs in its entirety and writes the output results to file.
- In the `test` case, the process runs in its entirety and compares the results to the output already existing in the designated output file.

Because the computational cost is extensive for a full run of the process of assembling TMCs, the $\mathbf{(5 \rightarrow 6)}$ process step applies a `demo` value for the argument `--xtent`, merely demonstrating the TMC assembly process step with a small selection of ligands for the monodentate or bidentate case.

## The Bash Files

The Bash files are to be run from the root of the repository, with the Conda environment `tmcinvdes` active. See [environment](environment/README.md).

For example:

```
sh bash/screen_generated_ligands-bidentate.sh
```

This will recapitulate the screening of generated ligands from `uncond_bi` and verify that the results are reproducible.

## Workflow Structure

With reference to the numbered stages $\mathbf{x}$, $\mathbf{y}$ described in [datasets/README.md](datasets/README.md), in the following we list the script filenames associated with for each process script $\mathbf{(x \rightarrow y)}$.

For example, the first bullet item below is labeled $\mathbf{(0 \rightarrow 1)}$ to designate that the process script takes dataset 0 (i.e., [tmQMg-L](https://github.com/hkneiding/tmQMg-L/)) as input and creates as output the datasets in stage 1 described in [datasets/README.md](datasets/README.md).

Any currently outstanding scripts will be added soon.

- $\mathbf{(0 \rightarrow 1)}$: Select initial ligands and get their encoded SMILES representation.
  - `tmcinvdes/ligand_generation/get_encoded_smiles.py`
  - `bash/get_encoded_smiles-monodentate.sh`
  - `bash/get_encoded_smiles-bidentate.sh`
- $\mathbf{(1 \rightarrow 2)}$: Train unconditional JT-VAE models. This is done with code in the separate JT-VAE repository [JT-VAE-tmcinvdes](https://github.com/Strandgaard96/JT-VAE-tmcinvdes/) applied to the `.txt` files in $1$.
- $\mathbf{(2 \rightarrow 3)}$: Unconditionally generate ligands with trained unconditional JT-VAE models. See our [JT-VAE-tmcinvdes](https://github.com/Strandgaard96/JT-VAE-tmcinvdes/) repository.
- $\mathbf{(3 \rightarrow 4)}$: Process and screen unconditionally generated ligands to extract only novel ligands.
  - `tmcinvdes/ligand_generation/screen_generated_ligands.py`
  - `bash/screen_generated_ligands-monodentate.sh`
  - `bash/screen_generated_ligands-bidentate.sh`
- $\mathbf{(4 \rightarrow 5)}$: Select $X\mathrm{k}$  smallest novel ligands.
  - `tmcinvdes/ligand_generation/select_Xk_smallest_ligands.py`
  - `bash/select_Xk_smallest_ligands-monodentate.sh`
  - `bash/select_Xk_smallest_ligands-bidentate.sh`
- $\mathbf{(5 \rightarrow 6)}$: Assemble TMC structures as XYZ blocks in concatenated output file.
  - `tmcinvdes/structure_generation/assemble_tmcs.py`
  - `bash/assemble_tmcs-monodentate-uncond_mono-min15k.sh`
  - `bash/assemble_tmcs-bidentate-uncond_bi-min10k.sh`
- $\mathbf{(6 \rightarrow 7)}$: Label TMCs by ORCA calculations on TMCs. The ORCA calculations are done on HPC, the present Python file only parses the ORCA output files.
  - `tmcinvdes/quantum_chemistry/ORCA/parse_orca_to_labels.py`
  - `bash/parse_orca_to_labels-monodentate-uncond_mono-min15k.sh`
  - `bash/parse_orca_to_labels-bidentate-uncond_bi-min10k.sh`
- $\mathbf{(7 \rightarrow 8a \text{ and } 8b)}$: Exclude outliers and include the rest of labeled TMCs.
  - `tmcinvdes/analysis/exclude_outliers.py`
  - `bash/exclude_outliers-monodentate.sh`
  - `bash/exclude_outliers-bidentate.sh`
- $\mathbf{(8b \rightarrow 9)}$: Train conditional JT-VAE models with ligands labeled by their TMCs' DFT properties. See our [JT-VAE-tmcinvdes](https://github.com/Strandgaard96/JT-VAE-tmcinvdes/) repository.
- $\mathbf{(8b \rightarrow 10)}$: Sample ligands for optimization.
  - `tmcinvdes/ligand_generation/sample_ligands_to_optimize.py`
  - `bash/sample_ligands_to_optimize-monodentate.sh`
- $\mathbf{(10 \rightarrow 11)}$: Optimize sampled ligands. See our [JT-VAE-tmcinvdes](https://github.com/Strandgaard96/JT-VAE-tmcinvdes/) repository.
- $\mathbf{(11 \rightarrow 12)}$: Assemble TMC structures as XYZ blocks in concatenated output file. (New Bash files but re-used Python script.)
  - `tmcinvdes/structure_generation/assemble_tmcs.py` *(Re-used from $\mathbf{(5 \rightarrow 6)}$.)*
  - `bash/assemble_tmcs-monodentate-cond_mono-sampled_optimized.sh`
- $\mathbf{(12 \rightarrow 13)}$: Label TMCs by ORCA calculations on TMCs. The ORCA calculations are done on HPC, the present Python file only parses the ORCA output files. (New Bash files but re-used Python script.)
  - `tmcinvdes/quantum_chemistry/orca/parse_orca_to_labels.py` *(Re-used from * $\mathbf{(6 \rightarrow 7)}$ *.)*
  - `bash/parse_orca_to_labels-monodentate-cond_mono-sampled_optimized.sh`

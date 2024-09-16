# Datasets and Files

The datasets are organized to be legible and consistent across the parallel workflows for monodentates and bidentates.

## Naming Conventions

The dataset subdirectories represent different fixed **stages** of the workflow, where a source stage is mapped to a target stage by some scripted **process step**.

- Ensure parallel, logical and intuitive naming of each stage.
- Append a hyphen and a term to the filename prefix for each downstream minor stage, so the filename prefix represents a starting point and a hyphen-separated sequence of processes applied.
- Major stages condense the process history of the previous major stage and its subsequent minor stages into a new, shorter name.

## Stages

Throughout we distinguish the parallel workflow stages with `_mono` for monodentate ligands and `_bi` for bidentate ligands.

Below we describe each stage in general terms, show the terms to be appended according to the naming conventions, list the filenames, and for `.csv` files we list the key column names.

Note that the 13 major bullet points below correspond to the directories into which we have prepared the data files. See the detailed description of the process steps between the stages and associated scripts [here](../DETAILS.md).

- **1. Select ligand subsets:** `_mono` or `_bi`
  - Ligands are selected from the zeroth stage [tmQMg-L](https://github.com/hkneiding/tmQMg-L) and are used to obtain encoded SMILES representations for use in the JT-VAE models.
  - Filenames:
    - `tmQMg-L_mono.csv`
    - `tmQMg-L_bi.csv`
  - Columns:
    - `tmQMg-L ligand ID`: str of the form `ligand-N-0` where `N` is an integer.
    - `Decoded SMILES`: str
    - `Connection IDs`: list\[int\]
    - `Encoded SMILES`: str
  - Usage: We keep a parallel (encoded SMILES only, no header) `.txt` file since that is what the current JT-VAE code expects.
- **2. Train unconditional generative models:**
  - Trained model files:
    - `uncond_mono.iter-10000`
    - `uncond_bi.iter-10000`
- **3. Unconditionally generate ligands:** `-raw50k`
  - Ligands generated using the trained unconditional JT-VAE with model as indicated by the respective beginning prefix in the filename.
  - Filenames:
    - `uncond_mono-raw50k.txt`
    - `uncond_bi-raw50k.txt`
  - The file format is left as `.txt` since this is the format produced by the current JT-VAE code. One encoded SMILES per row, across 50 000 rows.
- **4. Screen generated ligands:** `-novel`
  - Filenames:
    - `uncond_mono-raw50k-novel.csv`
    - `uncond_bi-raw50k-novel.csv`
  - Columns:
    - `Heavy atom count` (or `Atom count`) if monodentate (or bidentate), respectively: int
    - `Encoded SMILES`: str
    - `Decoded SMILES`: str
    - `Connection IDs`: list\[int\]
    - `Coordination environment`: list\[str\]
- **5. Rank by size, select smaller subset:** `-min15k` or `-min10k`
  - Note: we make a contraction here for brevity, taking `-minXk` to imply `-raw50k-novel-minXk`.
  - While the previous stage did not have an inherent ordering of ligands, this stage does, and we introduce identities for the ligands.
  - Difference: monodentate ligands are ordered by heavy atoms while bidentate ligands are ordered by total number of atoms, including hydrogen.
  - Filenames:
    - `uncond_mono-min15k.csv`
    - `uncond_bi-min10k.csv`
  - Columns:
    - `Ligand ID`: str of the form `uncond_bi-min15k-N` or `uncond_bi-min10k-N` where N is the integer of the ligand within this ligand set.
    - `Heavy atom count` (or `Atom count`) if monodentate (or bidentate), respectively: int
    - `Encoded SMILES`: str
    - `Decoded SMILES`: str
    - `Connection IDs`: list\[int\]
    - `Coordination environment`: list\[str\]
- **6. Assemble TMCs as XYZ structure files:** `-TMC`
  - Only the successfully generated XYZ files are used for DFT calculations in the next step, producing stage 7.
  - For `_bi`, in the XYZ block comment line we include the ligand ID and further distinguish using `-cis` and `-trans` based on structure generation input, i.e., intended isomerism.
  - For these TMCs, all the XYZ blocks are concatenated into one file for `mono` and one for `bi`, respectively.
  - Filenames:
    - `uncond_mono-min15k-TMC.xyz`
    - `uncond_bi-min10k-TMC.xyz`
- **7. Label ligands with DFT calculations on TMCs geometry optimizations, HOMO-LUMO gap, and metal center charge:** `-labeled`
  - The XYZ files generated from the ligand subsets are `labeled` with DFT properties calculated using ORCA, and the results were associated with the respective ligands.
    - For bidentate ligand-based TMCs, DFT calculations are performed only on those homoleptic TMCs from stage 6 where both isomer forms were successfully assembled for a given ligand.
  - We make another unambiguous contraction for brevity, taking `-labeled` to mean `-TMC-labeled`.
  - Filenames:
    - `uncond_mono-min15k-labeled.csv`
    - `uncond_bi-min10k-labeled.csv`
  - Columns:
    - `Ligand ID`: str of the form `uncond_bi-min15k-N` or `uncond_bi-min10k-N` where `N` is the integer of the ligand within this ligand set, as established in the `-minXk` step.
    - `Isomer` (bidentate only): str, `"cis"` or `"trans"` for bidentate ligands.
    - `Metal center charge`: float
    - `HOMO-LUMO gap (Eh)`: float
    - `XYZ`: str, a complete XYZ block for the TMC in question, including line breaks `\n` that are automatically handled (escaped/unsecaped) by Pandas when reading/writing from/to file.
- **8. Exclude (8a.) outliers with IsolationForest with respect to HOMO-LUMO gap and metal center charge, and include (8b.) remaining ligands:** `-excluded` or `-included`
  - Difference: ligands were excluded by a `contamination` argument of 1.0% for monodentates and 1.5% for bidentates.
  - Filenames, CSV files:
    - `uncond_mono-min15k-labeled-included.csv`
    - `uncond_bi-min10k-labeled-included.csv`
    - `uncond_mono-min15k-labeled-excluded.csv`
    - `uncond_bi-min10k-labeled-excluded.csv`
  - CSV columns: exactly as previous stage, i.e., as `...-labeled.csv` files.
  - Filenames, XYZ files:
    - `uncond_mono-min15k-labeled-excluded.xyz`
    - `uncond_bi-min10k-labeled-excluded.xyz`
    - Note that we keep the manual annotations made in the comment line of each XYZ block among the excluded outlier TMCs.
- **9. Train conditional generative models with labeled and included ligand data:** `cond`
  - Trained model files:
    - `cond_mono.epoch-149`
- **10. Sample ligands for optimization from the conditional model's own training set:** `-sampled_for_{conditional model name}`
  - These ligands are sampled to be used as starting points for the production of conditionally generated ligands.
  - Filenames:
    - `uncond_mono-min15k-labeled-included-sampled_for_cond_mono.csv`
  - Columns:
    - `Ligand ID`: str of the form `uncond_bi-min15k-N` or `uncond_bi-min10k-N` where `N` is the integer of the ligand within this ligand set as established in the `-minXk` step.
    - `Encoded SMILES`: str, retrieve from `-minXk` step.
    - `Metal center charge`: float
    - `HOMO-LUMO gap (Eh)`: float
    - `XYZ`: str, a complete XYZ block for the TMC in question.
- **11. Optimize sampled ligands:** `{model_name}-sampled_optimized`
  - At this point, the results are conditionally generated, so filenames can begin with the conditional model name.
  - Filename:
    - `cond_mono-sampled_optimized.csv`
  - Columns:
    - `Label ID`: int, from stage 7.
    - `Original ligand ID`: str
    - `Optimized ligand ID`: str, of the form `{model_name}-{optimization_direction}_optimized-{original_ligand_id}` where `optimization_direction` is `O1` through `O2` based on the compass directions in Table 2 of the preprint.
    - `Original decoded SMILES`: str
    - `Sampling region`: str, from which the original ligand was sampled.
      - Values:
        - `"Center"`
        - `"High metal center charge"`
        - `"Low HOMO-LUMO gap (eV)"`
        - `"High HOMO-LUMO gap (eV)"`
        - `"Low metal center charge"`
    - `Optimization objective`: a simplified expression corresponding to the objective function used to optimize the ligand. Here $\epsilon$ is the HOMO-LUMO gap (eV) and $q_\text{Ir}$ is the metal center charge.
      - Values:
        - $\epsilon$
        - $q_\text{Ir}$
        - $\epsilon + q_\text{Ir}$
        - $\epsilon - q_\text{Ir}$
    - `Minimization`: bool, whether to minimize or maximize the optimization objective.
    - `Original encoded SMILES`: str
    - `Optimized encoded SMILES`: str
    - `Tanimoto smiliarity`: float, comparison of original and optimized ligand.
    - `Original metal center charge`: float from stage 7.
    - `Original HOMO-LUMO gap (Eh)`: float from stage 7.
    - `Original HOMO-LUMO gap (eV)`: float from stage 7.
- **12. Generate XYZ structures of homoleptic TMCs using molSimplify:** `-TMC` (re-used)
  - Filename:
    - `cond_mono-sampled_optimized-TMC.xyz`
- **13. Label ligands with DFT calculations on TMCs geometry optimizations, HOMO-LUMO gap, and metal center charge:** `-labeled` (re-used)
  - We make the same unambiguous contraction for brevity, taking `-labeled` to mean `-TMC-labeled`.
  - Filename:
    - `cond_mono-sampled_optimized-labeled.csv`
  - Columns: as the `-sampled_optimized.csv` files above at stage 11, except:
    - `Optimized metal center charge`: float
    - `Optimized HOMO-LUMO gap (eV)`: float
    - `Original XYZ`: str, from stage 7.
    - `Optimized XYZ`: str

## Pandas DataFrames/CSV Conventions

At many stages in the workflow, we read from and write to .csv files, while the Python functions use Pandas to process data in DataFrames. This is a note on our normative conventions:

Here is how we write a DataFrame to file:

```Python
import pandas as pd

df = pd.DataFrame(...)

filepath = "foo"
df.to_csv(filepath, header=True, index=False)
```

This assumes our code explicitly handles column name ordering and internal IDs for different instances on rows, to be able to safely discard the automatic Pandas index information.

We can then later simply read from file:

```Python
filepath = "foo"
df = pd.read_csv(filepath)
```

# Datasets and Files

The datasets are organized to be legible and consistent across the parallel workflows for monodentates and bidentates.


## Naming Conventions

The dataset files represent different fixed **stages** of the workflow, where a source stage is piped to a target stage by some script process (function).

 - Ensure parallel, logical and intuitive naming of each stage.
 - Append a hyphen and a term to the filename prefix for each downstream minor stage, so the filename represents a starting point and a hyphen-separated sequence of processes applied.
 - Major stages condense the process history of the previous major stage and its subsequent minor stages into a new, shorter name.


## Stages

Throughout we distinguish the parallel workflow stages with `_mono` for monodentate ligands and `_bi` for bidentate ligands.
Note that the 14 major bullet points below correspond to the directories into which we should prepare the data files.

 - **1. Select ligand subsets:** `_mono` or `_bi`
   - [tmQMg-L](https://github.com/hkneiding/tmQMg-L) is the starting point, using `ligands_descriptors.csv`,  `ligands_fingerprints.csv`, `ligands_misc_info.csv`.
   - Script: these data files are to be made reproducible/obtainable by `tmcinvdes/structure_generation/select_ligands.py` via `bash/select-tmQMg-L-monodentates.sh` / `bash/select-tmQMg-L-bidentates.sh`.
   - Filenames:
      - `tmQMg-L_mono.csv`
      - `tmQMg-L_bi.csv`
   - Columns:
       - `tmQMg-L ligand ID`: str of the form `ligand-N-0` where `N` is an integer.
       - `Canonical SMILES`: str
       - `Connection IDs`: int, list[int,[list[int]]]
       - `Enriched SMILES`: str
   - Usage: We keep a parallel (enriched SMILES only, no header) `.txt` file since that is what the current JT-VAE code expects, and we use the same concatenated prefix (+`-vocab`) for the filename of the vocab text file under `vocabs/` (e.g., `vocabs/tmQMg-L_mono-vocab.txt`) in the JT-VAE repository.
 - **2. Train unconditional generative models:** 
     - Trained model files:
         - `uncond_mono.iter-10000` 
         - `uncond_bi.iter-10000` 
 - **3. Unconditionally generate ligands:** `-raw50k`
     - Ligands generated using the trained unconditional JT-VAE with model as indicated by the respective beginning prefix in the filename.
     - Filenames:
         - `uncond_mono-raw50k.txt`
         - `uncond_bi-raw50k.txt`
     - The file format is left as `.txt` since this is the format produced by the current JT-VAE code. One enriched SMILES per row, across 50 000 rows.
 - **4. Deduplicate, filter for novelty, process ligands, and screen out incorrect coordination modes:** `-novel`
     - Filenames:
         - `uncond_mono-raw50k-novel.csv`
         - `uncond_bi-raw50k-novel.csv`
     - Columns:
         - `Heavy atom count` (or `Atom count`) if monodentate (or bidentate), respectively: int
         - `Enriched SMILES`: str
         - `Canonical SMILES`: str
         - `Connection IDs`: int, list[int,[list[int]]]
 - **5. Rank by size and select smaller subset:** `-min15k` or `-min10k`
     - Note: we make a contraction here for brevity, taking `-minXk` to imply `-raw50k-novel-minXk`.
     - Note: while the previous stage did not have an inherent ordering of ligands, this stage does, and we introduce identities for the ligands.
     - Difference: monodentate ligands are ordered by heavy atoms while bidentate ligands are ordered by total number of atoms, including Hydrogen.
     - Filenames:
         - `uncond_mono-min15k.csv`
         - `uncond_bi-min10k.csv`
     - Columns:
         - `Ligand ID`: str of the form `uncond_bi-min15k-N` or `uncond_bi-min10k-N` where N is the integer of the ligand within this ligand set.
         - `Heavy atom count` (or `Atom count`) if monodentate (or bidentate), respectively: int
         - `Enriched SMILES`: str
         - `Canonical SMILES`: str
         - `Connection IDs`: int, list[int,[list[int]]]
 - **6. Generate XYZ structures of homoleptic TMCs using molSimplify:** `-TMC`
     - Only the successfully molSimplify-generated XYZ files are used for DFT calculations in the next step, producing stage 7.
     - For `_bi`, in the XYZ block comment line we include the ligand ID and further distinguish using `-cis` and `-trans` based on structure generation input, i.e., intended isomerism.
     - For these TMCs, all the XYZ blocks are concatenated into one file for `mono` and one for `bi`, respectively.
     - Filenames:
         - `uncond_mono-min15k-TMC.xyz`
         - `uncond_bi-min10k-TMC.xyz`
 - **7. Label ligands with DFT calculations on TMCs geometry optimizations, HOMO-LUMO gap, and metal center charge:** `-labeled`
     - The XYZ files generated from the ligand subsets are `labeled` with DFT properties calculated using ORCA, and the results were associated with the respective ligands.
         - For bidentate ligand-based TMCs, DFT calculations are performed only on the subset from stage 6 where a ligand is present in both isomer forms om the homoleptic TMCs. 
     - We make another unambiguous contraction for brevity, taking `-labeled` to mean `-TMC-labeled`.
     - Filenames:
         - `uncond_mono-min15k-labeled.csv`
         - `uncond_bi-min10k-labeled.csv`
     - Columns:
         - `Ligand ID`: str of the form `uncond_bi-min15k-N` or `uncond_bi-min10k-N` where `N` is the integer of the ligand within this ligand set, as established in the `-minXk` step.
         - `Isomer`: str, either empty string `""` for monodentate ligands, or else `"cis"` or `"trans"` for bidentate ligands.
         - `Metal center charge`: float
         - `HOMO-LUMO gap (Eh)`: float
         - `XYZ`: str, a complete XYZ block for the TMC in question, including line breaks `\n` that are automatically handled (escaped/unsecaped) by Pandas when reading/writing from/to file.
         -    - NB: make sure the comment line of the XYZ block contains the correct ligand ID (plus `-cis` or `-trans` for bidentate ligands).
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
         - Note that we keep the manual annotations made in the comment line of each XYZ block.
 - **9. Train conditional generative models with labeled and included ligand data:** remap names of models to `cond_mono`, `cond_bi_cis`, `cond_bi_trans`(, and eventually `cond_bi_cistrans` for the combined isomers model with symmetry check).
     - Trained model files:
         - `cond_mono.epoch-149`
         - *(`cond_bi_cis.epoch-149`)*
 - **10. Sample ligands for optimization from the conditional model's own training set:** `-sampled_for_{conditional model name}`
     - These ligands are sampled to be used as starting points for the production of conditionally generated ligands.
     - Filenames:
         - `uncond_mono-min15k-labeled-included-sampled_for_cond_mono.csv`
         - *(`uncond_bi-min10k-labeled-included-sampled_for_cond_bi_cis.csv`)*
     - Columns:
         - `Ligand ID`: str of the form `uncond_bi-min15k-N` or `uncond_bi-min10k-N` where `N` is the integer of the ligand within this ligand set as established in the `-minXk` step.
         - `Enriched SMILES`: str, retrieve from `-minXk` step.
         - `Isomer`: str, either empty string `""` for monodentate ligands, or else `"cis"` or `"trans"` for bidentate ligands.
         - `Metal center charge`: float
         - `HOMO-LUMO gap (Eh)`: float
         - `XYZ`: str, a complete XYZ block for the TMC in question.
 - **11. Optimize sampled ligands:** `{model name}-sampled_optimized`
     - At this point, the results are conditionally generated, so filenames can begin with the conditional model name.
     - Filenames:
         - `cond_mono-sampled_optimized.csv`
         - *(`cond_bi_cis-sampled_optimized.csv`)*
     - Columns:
         - TBD
 - **12. Generate XYZ structures of homoleptic TMCs using molSimplify:** `-TMC` (re-used)
     - Filenames:
         - `cond_mono-sampled_optimized-TMC.xyz`
         - *(`cond_bi_cis-sampled_optimized-TMC.xyz`)*
 - **13. Label ligands with DFT calculations on TMCs geometry optimizations, HOMO-LUMO gap, and metal center charge:** `-labeled` (re-used)
     - Filenames:
         - `cond_mono-sampled_optimized-TMC-labeled.csv`
         - *(`cond_bi_cis-sampled_optimized-TMC-labeled.csv`)*
    - Columns: same as the `-labeled.csv` files above at stage 7.


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

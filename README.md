# tmcinvdes

Repository for the preprint ["Deep Generative Model for the Dual-Objective Inverse Design of Metal Complexes."](https://doi.org/10.26434/chemrxiv-2024-mzs7b)

`tmcinvdes` will provide code and data from our research around Junction Tree Variational Autoencoder (JT-VAE) models that can generate metal ligands for transition metal complexes (TMCs). TMCs are used in industrial catalytic processes, anticancer therapies, and energy transformations. By labeling ligands with DFT-calculated target properties, conditional generative models can be trained and harnessed, optimizing metal ligands directionally in the target property space to discover novel **TMC**s by **inv**erse **des**ign.

<img align="center" src="concept_overview.png" alt="Inverse Design of Metal Complexes" width="800"/>

## Code

**The code used to train the JT-VAE models and generate ligands are found at: [JT-VAE-tmcinvdes](https://github.com/Strandgaard96/JT-VAE-tmcinvdes).**

Follow the environment [setup instructions](/environment/README.md) to install all the dependencies of the code in this repository.
The code relies on a local download of the [tmQMg-L](https://github.com/hkneiding/tmQMg-L.git) repository and must therefore be cloned into a local folder.

[Ligand generation](/tmcinvdes/ligand_generation)

Contains the code used to create the JT-VAE training sets.

[Structure generation](/tmcinvdes/structure_generation)

Contains the code to assemble TMCs from ligands.

[Quantum chemistry](/tmcinvdes/quantum_chemistry)

Contains the ORCA input file and parser scripts to label the generated ligands with the DFT-calculated properties of their homoleptic TMCs.

[Analysis](/tmcinvdes/analysis)

Contains the SA score code and the code for excluding outliers.

For a description of the workflow structure see the [detailed workflow](DETAILS.md).

### Citation

If you find our work useful, please cite our article:

```
@article{Strandgaard:2024:ChemRxiv,
  author    = {Magnus Strandgaard, Trond Linjordet, Hannes Kneiding, Arron Burnage, Ainara Nova, Jan Halborg Jensen, David Balcells}
  title     = {Deep Generative Model for the Dual-Objective Inverse Design of Metal Complexes},
  journal   = {ChemRxiv},
  volume    = {preprint},
  doi       = {10.26434/chemrxiv-2024-mzs7b-v2},
  year      = {2024},
}
```

### Work in Progress

This repo will be gradually updated with code and data as the preprint goes through review.

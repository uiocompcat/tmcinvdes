# tmcinvdes

Public repository for the preprint ["Deep Generative Model for the Dual-Objective Inverse Design of Metal Complexes."](https://doi.org/10.26434/chemrxiv-2024-mzs7b)

`tmcinvdes` will provide code and data from our research around junction tree variational autoencoder (JT-VAE) models that can generate metal ligands for transition metal complexes (TMCs). TMCs are used in industrial catalytic processes, anticancer therapies, and energy transformations. By labeling ligands with DFT-calculated target properties, conditional generative models can be trained and harnessed, optimizing metal ligands directionally in the target property space to discover novel, useful **TMC**s by **inv**erse **des**ign.

<img align="center" src="concept_overview.png" alt="Inverse Design of Metal Complexes" width="800"/>

The code used to train the JT-VAE models and generate ligands can be found at: [FastJTNNpy3](https://github.com/Strandgaard96/FastJTNNpy3).

## Code

Follow the environment [setup instructions](/environment/README.md) to install all the dependencies of the code in this repository.
The code also relies on a local download of the [tmQMg-L](https://github.com/hkneiding/tmQMg-L.git) repository.

[Ligand generation](/tmcinvdes/ligand_generation)

Contains the code used to create the JT-VAE training sets.

[Structure generation](/tmcinvdes/structure_generation)

Contains the code to assemble TMCs from ligands.

[Quantum chemistry](/tmcinvdes/quantum_chemistry)

Contains the ORCA input file and parser scripts to label the generated ligands with the DFT-calculated properties of their homoleptic TMCs.

[Analysis](/tmcinvdes/analysis)

Contains the code for analyzing and processing data.


For nominal runs of the present code, shared both as explicit usage examples and to test reproducibility, see the [detailed workflow](DETAILS.md).

### Work in Progress

This repository represents work in progress, meaning that data and code will be gradually added.

# tmcinvdes

Public repository for the preprint ["Deep Generative Model for the Dual-Objective Inverse Design of Metal Complexes."](https://doi.org/10.26434/chemrxiv-2024-mzs7b)

`tmcinvdes` will provide code and data from our research around junction tree variational autoencoder (JT-VAE) models that can generate metal ligands for transition metal complexes (TMCs). TMCs are used in industrial catalytic processes, anticancer therapies, and energy transformations. By labeling ligands with DFT-calculated target properties, conditional generative models can be trained and harnessed, optimizing metal ligands directionally in the target property space to discover novel, useful **TMC**s by **inv**erse **des**ign.

<img align="center" src="concept_overview.png" alt="Inverse Design of Metal Complexes" width="800"/>

The fork used to train the models can be found at: [FastJTNNpy3](https://github.com/Strandgaard96/FastJTNNpy3)

## Work in Progress

This repository represents work in progress, meaning that data and code will be gradually added.

One great way to stay informed is to make sure you are logged in to GitHub and then click the `Watch` button with the eye-icon, and choose appropriate options.

For example, if only major updates are of interest: click `Watch` -> `Custom` -> check the box `Releases` -> `Apply`.

## Datasets

The first priority is to release datasets associated with the preprint. These will be organized into numbered subdirectories under the directory `datasets`. The dataset subdirectories are numbered in sequence to clarify the workflow. Note that some files and subdirectories described in `datasets/README.md` will be added in the future.

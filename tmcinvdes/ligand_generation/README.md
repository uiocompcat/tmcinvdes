# Ligand generation

This folder contains code for generating the training data for the JT-VAE

To generate the training sets run the following script. Datafiles from the tmQMg-L need to be in a local folder and supplied as input to the -i argument.

```
python select_initial_ligands.py -i ~/git/tmQMg-L -o output_dir
```

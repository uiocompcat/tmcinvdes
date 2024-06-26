# Ligand Generation

This folder contains code for selecting the training data for the unconditional JT-VAE models.

To generate the training sets run the following script. The [tmQMg-L](https://github.com/hkneiding/tmQMg-L) repository needs to be in a local folder and supplied as input to the -i argument.

```
python get_substituted_smiles.py -i ~/git/tmQMg-L -o output_dir -d monodentate
```

## RDKIT version

Be aware that minor differences in RDKit version will lead to different outputs of `get_substituted_smiles.py`.

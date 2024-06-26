# How to set up conda environment

## 1. Create Conda Environment

Make sure you have conda installed and from this directory run this command:

```
conda env create -f tmcinvdes.yml
```

Note that this installs the package [uxtbpy](https://github.com/hkneiding/uxtbpy). This is suitable for running on shared nodes where you may not have the rights to add symbolic links arbitrarily or should not because it could affect others' workflows.

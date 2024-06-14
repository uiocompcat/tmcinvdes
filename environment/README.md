# How to Set Up

This document describes the setup of the expected environment.

The [first section](/environment/README.md#1-install-pre-requisite-binaries) describes how to install the pre-requisite binaries using (A) the [express approach](/environment/README.md#express-installation) using a Bash, or (B) how to do the same steps [manually](/environment/README.md#step-wise-installation) (documented for troubleshooting).

The [second section](/environment/README.md#2-create-conda-environment) describes how to create the conda envionment.

It is recommended to install the binaries in Section 1 using only the express approach, then install the conda environment as described in Section 2.


## 1. Install Pre-Requisite Binaries

### Express Installation

Starting from the project root, navigate to the `environment/` directory and run the installation script `install_xtb_stda.sh` as follows:
```
cd environment
bash install_xtb_stda.sh
```
This script download and make available all binaries required for running full xTB-sTDA calculations. Namely, the individual components are:

 - xtb (The base xtb program used to optimize structures and calculate quantum properties.)
 - xtb4stda (xtb-stda adapter program used to obtain a wavefunction for a given molecule to be used as input to the stda program.)
 - stda (stda program used to calculate excited states according to the formulation of the simplified Tamm-Dancoff Approximation.)
 - g_spec (Auxilliary program used to determine line and broaden spectra based on stda output.)

Note that this will create a subdirectory `xtb_stda_binaries/` in your current working directory to store the binaries. If you want this directory to be located somewhere else than the `environments/` directory you should first navigate to your target directory and then run the installation script with:
```
bash /absolute/path/to/install_xtb_stda.sh
```

### Step-wise Installation

The following is an alterative to the express installation of the pre-requisite binaries.
Create and navigate into a directory to store the binaries associated to this project.

#### Install xTB (for Linux x86)
Next, do this to install xTB:
```
wget https://github.com/grimme-lab/xtb/releases/download/v6.4.1/xtb-6.4.1-linux-x86_64.tar.xz
tar -xvf xtb-6.4.1-linux-x86_64.tar.xz
rm xtb-6.4.1-linux-x86_64.tar.xz
ln -s xtb-6.4.1/bin/xtb xtb
chmod +x xtb
```

#### Install STDA (for Linux x86)

In the same directory, do this to install sTDA:

```
wget https://github.com/grimme-lab/stda/releases/download/v1.6.3/xtb4stda
wget https://github.com/grimme-lab/stda/releases/download/v1.6.3/stda_v1.6.3
wget https://github.com/grimme-lab/stda/releases/download/v1.6.3/g_spec

sudo chmod +x xtb4stda
sudo chmod +x stda_v1.6.3
sudo chmod +x g_spec
```

#### Make the binaries available on the whole system
In the same directory run:
```
echo export PATH=$PATH:$(pwd) >> ~/.bashrc
source ~/.bashrc
```

#### Download the parameter files
Finally, download the parameter files into the user root:
```
wget https://raw.githubusercontent.com/grimme-lab/xtb4stda/master/.param_stda1.xtb -P ~/
wget https://raw.githubusercontent.com/grimme-lab/xtb4stda/master/.param_stda2.xtb -P ~/
```

## 2. Create Conda Environment

Once the pre-requisite binaries are installed as per Section 1, ensure that you have an active conda (base) environment, then in the `environment/` directory run this command:
```
conda env create -f tmcinvdes.yml
```
Note that this installs the package [uxtbpy](https://github.com/hkneiding/uxtbpy) which depends on the xTB installed in the first step. This is suitable for running on shared nodes where you may not have the rights to add symbolic links arbitrarily or should not because it could affect others' workflows.

Instantiate a local molSimplify ligand library, assuming you do not already have one:
```
mkdir molsimp_local
cd molsimp_local
yes "$(pwd)" | molsimplify -ligadd "S=CC(=S)Br" -ligname uncond_bi-min10k-1 -ligcon 1,4 -skipANN True
```
You can create the directory for the molSimplify ligand library anywhere on your machine.

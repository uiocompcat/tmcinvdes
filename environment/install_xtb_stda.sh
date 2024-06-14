#!/bin/bash -l

# download binaries
wget https://github.com/grimme-lab/xtb/releases/download/v6.4.1/xtb-6.4.1-linux-x86_64.tar.xz
wget https://github.com/grimme-lab/stda/releases/download/v1.6.3/xtb4stda
wget https://github.com/grimme-lab/stda/releases/download/v1.6.3/stda_v1.6.3
wget https://github.com/grimme-lab/stda/releases/download/v1.6.3/g_spec

# extract xtb
tar -xvf xtb-6.4.1-linux-x86_64.tar.xz
rm xtb-6.4.1-linux-x86_64.tar.xz

# symbolic link for xtb
ln -s xtb-6.4.1/bin/xtb xtb

# add execution permissions
chmod +x xtb
chmod +x xtb4stda
chmod +x stda_v1.6.3
chmod +x g_spec

# download parameter files
wget https://raw.githubusercontent.com/grimme-lab/xtb4stda/master/.param_stda1.xtb -P ~/
wget https://raw.githubusercontent.com/grimme-lab/xtb4stda/master/.param_stda2.xtb -P ~/

# add to path variable
echo export PATH=$PATH:$(pwd) >> ~/.bashrc
source ~/.bashrc

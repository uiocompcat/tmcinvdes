#!/bin/bash
python -m tmcinvdes.structure_generation.assemble_tmcs -d bidentate \
                                                       -i datasets/11_cond-optimized/cond_bi-isolated_ligands-sampled_optimized-to_verify.csv \
                                                       -o datasets/12_cond-TMC/cond_bi-isolated_ligands-sampled_optimized-TMC.xyz \
                                                       -p True \
                                                       -x full

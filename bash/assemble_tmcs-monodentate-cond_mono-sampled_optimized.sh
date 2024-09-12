#!/bin/bash
python -m tmcinvdes.structure_generation.assemble_tmcs -d monodentate \
                                                       -i datasets/11_cond-optimized/cond_mono-sampled_optimized.csv \
                                                       -o datasets/12_cond-TMC/cond_mono-sampled_optimized-TMC.xyz \
                                                       -p True \
                                                       -x demo

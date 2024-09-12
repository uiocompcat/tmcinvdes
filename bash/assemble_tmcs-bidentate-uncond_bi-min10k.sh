#!/bin/bash
python -m tmcinvdes.structure_generation.assemble_tmcs -d bidentate \
                                                       -i datasets/05_uncond_minXk/uncond_bi-min10k.csv \
                                                       -o datasets/06_uncond-TMC/uncond_bi-min10k-TMC.xyz \
                                                       -x demo

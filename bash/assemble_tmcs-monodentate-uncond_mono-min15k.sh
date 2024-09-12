#!/bin/bash
python -m tmcinvdes.structure_generation.assemble_tmcs -d monodentate \
                                                       -i datasets/05_uncond_minXk/uncond_mono-min15k.csv \
                                                       -o datasets/06_uncond-TMC/uncond_mono-min15k-TMC.xyz \
                                                       -x demo

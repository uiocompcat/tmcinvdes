#!/bin/bash
python -m tmcinvdes.quantum_chemistry.parse_orca_to_labels -d monodentate \
                                                           -i tmp_results/orca_out-cond_mono-optimized-TMC/ \
                                                           -o datasets/13_cond-labeled/cond_mono-sampled_optimized-labeled.csv \
                                                           -r datasets/11_cond-optimized/cond_mono-sampled_optimized.csv \
                                                           -t optimized \
                                                           -x test

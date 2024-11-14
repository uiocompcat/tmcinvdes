#!/bin/bash
python -m tmcinvdes.quantum_chemistry.parse_orca_to_labels -d bidentate \
                                                           -i tmp_results/orca_out-cond_bi-DFT-optimized-TMC/ \
                                                           -o datasets/13_cond-labeled/cond_bi-DFT-sampled_optimized-labeled.csv \
                                                           -r datasets/11_cond-optimized/cond_bi-DFT-sampled_optimized.csv \
                                                           -t optimized \
                                                           -x test

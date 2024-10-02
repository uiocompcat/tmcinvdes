#!/bin/bash
python -m tmcinvdes.quantum_chemistry.parse_orca_to_labels -d bidentate \
                                                           -i tmp_results/orca_out-uncond_bi-min10k-TMC/ \
                                                           -o datasets/07_uncond-labeled/uncond_bi-min10k-labeled.csv \
                                                           -s datasets/05_uncond-minXk/uncond_bi-min10k.csv \
                                                           -x test

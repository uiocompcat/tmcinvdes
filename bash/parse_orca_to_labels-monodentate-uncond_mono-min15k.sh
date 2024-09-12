#!/bin/bash
python -m tmcinvdes.quantum_chemistry.parse_orca_to_labels -d monodentate \
                                                           -i tmp_results/orca_out-uncond_mono-min15k-TMC/ \
                                                           -o datasets/07_uncond-labeled/uncond_mono-min15k-labeled.csv \
                                                           -x test

#!/bin/bash
python -m tmcinvdes.ligand_generation.sample_ligands_to_optimize -d bidentate \
                                                                 -i datasets/01_tmQMg-L-training_sets/ \
                                                                 -l isolated_ligands \
                                                                 -o datasets/10b_sampled_from_01/ \
                                                                 -s 01 \
                                                                 -x full

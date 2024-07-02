#!/bin/bash
python -m tmcinvdes.ligand_generation.screen_generated_ligands -d monodentate \
                                                               -i datasets/03_uncond-raw50k-generated_ligands \
                                                               -o datasets/04_uncond_novel \
                                                               -t datasets/01_tmQMg-L-training_sets \
                                                               -x test

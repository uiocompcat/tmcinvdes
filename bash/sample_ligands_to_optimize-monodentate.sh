#!/bin/bash
python -m tmcinvdes.ligand_generation.sample_ligands_to_optimize -d monodentate \
                                                                 -i datasets/08b_uncond-included/ \
                                                                 -o datasets/10_uncond-sampled_from_8b/ \
                                                                 -x test

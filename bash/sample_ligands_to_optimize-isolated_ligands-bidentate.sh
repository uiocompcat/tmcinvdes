#!/bin/bash
python -m tmcinvdes.ligand_generation.sample_ligands_to_optimize -d bidentate \
                                                                 -i datasets/08b_uncond-included/ \
                                                                 -l isolated_ligands \
                                                                 -o datasets/10_uncond-sampled_from_8b/ \
                                                                 -x test

#!/bin/bash
python -m tmcinvdes.ligand_generation.select_Xk_smallest_ligands -d bidentate \
                                                                 -i datasets/04_uncond_novel \
                                                                 -o datasets/05_uncond_minXk \
                                                                 -x test

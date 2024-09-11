#!/bin/bash
python -m tmcinvdes.ligand_generation.get_encoded_smiles -d bidentate -i ~/Documents/tmQMg-L \
                                                         -o datasets/01_tmQMg-L-training_sets \
                                                         -x test

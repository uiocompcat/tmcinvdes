#!/bin/bash
python -m tmcinvdes.structure_generation.prepare_xyz_for_dft \
                -i datasets/12_cond-TMC/cond_bi-isolated_ligands-sampled_optimized-TMC.xyz \
                -o tmp_results/orca_xyz-cond_bi-isolated_ligands-sampled_optimized-TMC

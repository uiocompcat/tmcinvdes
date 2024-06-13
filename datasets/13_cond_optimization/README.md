# Conditional optimization results

[conditional_optimization_mono_results.csv](conditional_optimization_mono_results.csv) contains the results from conditionally optimizing 160 monodentate ligands in 8 directions.

Columns:

- `property_to_optimize` denotes which property was optimized. 'first' is HOMO-LUMO gap and 'second' is Ir charge.
- `minimize` denotes whether the property was optimized or minimized.
- `starting_region` denotes where in DFT labeled space the prompt ligand was sampled from.

[279_novel_dft_verified](279_novel_dft_verified.csv) contains the results for 279 of the 1280 ligands that where successfully verified with DFT. The columns with the '\_new' suffix are the optimized values.

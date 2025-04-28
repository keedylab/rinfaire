# RINFAIRE example scripts

1. Generate a "multinetwork" containing individual networks for a set of structures:

run_generate_multi.sh

2. Generate a sum network for the entire multinetwork, including scaling and visualization:

run_analysis_sum.sh

3. Generate a sum network for a subset of the multinetwork, including scaling and visualization (alternative to step 2):

run_analysis_sum_subset.sh

4.Evaluate overlap between sum network and user-defined residues of interest, including statistical significance:

run_analysis_residues_of_interest.sh

Please note that additional example scripts in the directory 'work in progress' may not have been tested with the current code so should be used with caution. 

Input .pkl files for run_analysis_residues_of_interest.sh and run_analysis_sum_subset.sh can be accessed on Zenodo (refer to paper).

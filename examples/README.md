# RINFAIRE

1. Generate individual networks for a set of structures. Use a .txt file with directory information of each input file (refer to All_PTPs_input.txt in /examples).
run_generate_multi.sh

2. Subsetting can be performed using the flag --subset and specifying column name. A metadata file in csv format will be supplied in step 1 to achieve this.
run_analysis_sum_subset.sh

3. Use the MultiNetwork.pkl output to run sum analysis.
run_analysis_sum.sh

4. Additional analysis including # of residues in the network that are within 4A of some residues of interest can be performed. This data can be supplied as another csv file. Please output a pkl file in step 3 to serve as input here using the --output_pickle flag in step 3.
run_analysis_residues_of_interest

Additional examples are provided; however, please note that they have not been tested or implemented in the current code, so use them with caution.

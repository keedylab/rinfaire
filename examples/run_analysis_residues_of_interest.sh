inputdir='/your/input/dir/'
outputfilename='/your/output/dir/output_filename'
inputsetpath='/multirin/installation/dir/examples/ResiduesOfInterest_Fig6.csv

pipenv run python /multirin/installation/dir/analysis_residues_of_interest.py \
	${inputdir}SumNetwork.pkl \
	${outputfilename} \
	--input_set ${inputsetpath} \
	--col Sector_B \
	--find_significance reference.pdb \
	--n_iter_sig_test 5 \
	--cumulative_histogram \
	--no_normalize_by_total \


inputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/testruns/test23/'
outputfilename='/data/araju/ptpsinpdb/MultiRIN_Outputs/testruns/test23/SubNetwork_All'

pipenv run python ../scripts/generate_subset.py \
	${inputdir}MultiNetwork.pkl \
	'Mutation_Uniprot_Residue' \
	${outputfilename} \
	--make_discrete 50 0 300 \
	--group '0-50' \
	# --make_discrete 0.5 1.0 1.5 \
	# --group '1.0-1.5' \


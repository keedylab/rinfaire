
inputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run13/'
outputfilename='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run13/SubNetwork_All'

pipenv run python ../scripts/generate_subset.py \
	${inputdir}MultiNetwork.pkl \
	'Protein_Name' \
	${outputfilename} \
	# --group 'STEP' \


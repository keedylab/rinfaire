
inputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run13/'
outputfilename='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run13/SubNetwork_All'

pipenv run python ../scripts/generate_subset.py \
	${inputdir}MultiNetwork.pkl \
	'Resolution' \
	${outputfilename} \
	--make_discrete 0.5 1.0 1.5 \
	# --group '1.0-1.5' \


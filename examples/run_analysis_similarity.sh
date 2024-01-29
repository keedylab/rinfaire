
inputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run7/'
outputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run7/'

pipenv run python ../scripts/analysis_similarity.py \
	${inputdir}MultiNetwork.pkl \
	${outputdir} \
	--distance_metric deltacon


inputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run6/'
outputfilename='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run6/covarianceMatrix'

pipenv run python ../scripts/analysis_covariance.py \
	${inputdir}MultiNetwork.pkl \
	${outputfilename} \
	-r

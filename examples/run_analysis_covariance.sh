
inputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run7/'
outputfilename='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run7/covarianceMatrix'

pipenv run python ../scripts/analysis_covariance.py \
	${inputdir}MultiNetwork.pkl \
	${outputfilename} \
	-r \
	--remove_weak_edges 25


inputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run7/'
outputfilename='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run7/SumNetwork_V4'

pipenv run python ../scripts/analysis_sum.py \
	${inputdir}MultiNetwork.pkl \
	${outputfilename} \
	--remove_weak_edges 95

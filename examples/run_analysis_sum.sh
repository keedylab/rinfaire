
inputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run6/'
outputfilename='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run6/SumNetwork_V2'

pipenv run python ../scripts/analysis_sum.py \
	${inputdir}MultiNetwork.pkl \
	${outputfilename} \
	--remove_weak_edges 25

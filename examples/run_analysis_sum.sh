
inputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/testruns/test18/'
outputfilename='/data/araju/ptpsinpdb/MultiRIN_Outputs/testruns/test18/SumNetwork'

pipenv run python ../scripts/analysis_sum.py \
	${inputdir}MultiNetwork.pkl \
	${outputfilename}

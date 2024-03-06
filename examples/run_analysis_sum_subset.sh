
inputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run13/'
outputfilename='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run13/SumNetwork_Batch_V1'

pipenv run python ../scripts/analysis_sum.py \
	${inputdir}MultiNetwork.pkl \
	${outputfilename} \
    --subset 'Protein_Name' \
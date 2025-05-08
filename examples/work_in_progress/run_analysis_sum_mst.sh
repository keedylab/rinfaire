inputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run19/'
outputfilename='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run19/SumNetwork_MST'

pipenv run python ../scripts/analysis_sum.py \
	${inputdir}MultiNetwork.pkl \
	${outputfilename} \
	--seq_to_ref /data/araju/ptpsinpdb/FinalProcessing/Final/1SUG_qFit_chainA.pdb \
    --mst \
	--output_graph_info \
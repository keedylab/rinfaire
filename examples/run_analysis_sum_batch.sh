
inputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run20/'
outputfilename='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run20/SumNetwork_Removed95'

pipenv run python ../scripts/analysis_sum.py \
	${inputdir}MultiNetwork.pkl \
	${outputfilename} \
	--remove_weak_edges 95 \
	--seq_to_ref /data/araju/ptpsinpdb/FinalProcessing/Final/1SUG_qFit_chainA.pdb \
	--remove_subgraphs 5 \
	--detect_communities \
	--output_modularity \
	--n_communities 2 \
	--output_pickle 

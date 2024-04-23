
inputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run19/'
outputfilename='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run19/SumNetwork_NoRemoval'

pipenv run python ../scripts/analysis_sum.py \
	${inputdir}MultiNetwork.pkl \
	${outputfilename} \
	--remove_weak_edges 95 \
	--remove_subgraphs 5 \
	--seq_to_ref /data/araju/ptpsinpdb/FinalProcessing/Final/1SUG_qFit_chainA.pdb \
	--output_graph_info \
	# --detect_communities \
	# --output_modularity \
	# --n_communities 6

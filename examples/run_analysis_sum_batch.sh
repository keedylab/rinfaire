
inputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run39/'
outputfilename='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run39/SumNetwork_Struct_Scaling'

pipenv run python ../scripts/analysis_sum.py \
	${inputdir}MultiNetwork.pkl \
	${outputfilename} \
	--remove_weak_edges 0 \
	--scale_sum_network struct \
	# --seq_to_ref /data/araju/ptpsinpdb/FinalProcessing/Final/1SUG_qFit_chainA.pdb \
	# --remove_subgraphs 5 \
	# --detect_communities \
	# --output_modularity \
	# --n_communities 7 \
	# --output_pickle 

inputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run7/'
outputfilename='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run7/SumNetwork_Removed95_Cluster_Pickle'

pipenv run python ../scripts/analysis_sum.py \
	${inputdir}MultiNetwork.pkl \
	${outputfilename} \
	--remove_weak_edges 95 \
	-c \
	--n_communities 5 \
	--remove_subgraphs 10 \
	--output_modularity \
	--output_graph_info \
	--output_pickle \
	# --seq_to_ref /data/araju/ptpsinpdb/FinalProcessing/Final/3EAX_qFit_chainA.pdb \
	# -p

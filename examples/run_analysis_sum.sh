inputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run7/'
outputfilename='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run7/SumNetwork_Removed95_Cluster_V5'

pipenv run python ../scripts/analysis_sum.py \
	${inputdir}MultiNetwork.pkl \
	${outputfilename} \
	--remove_weak_edges 95 \
	-c \
	--n_communities 5 \
	--remove_subgraphs 10 \
	--output_modularity \
	--seq_to_ref /data/araju/ptpsinpdb/FinalProcessing/Final/6B8Z_qFit_chainA.pdb \
	# -p
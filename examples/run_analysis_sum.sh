inputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/testruns/test42/'
outputfilename='/data/araju/ptpsinpdb/MultiRIN_Outputs/testruns/test42/SumNetwork_Removed0_Info'

pipenv run python ../scripts/analysis_sum.py \
	${inputdir}MultiNetwork.pkl \
	${outputfilename} \
	--remove_weak_edges 0 \
	--scale_sum_network struct \
	--output_graph_info \
	# -c \
	# --n_communities 5 \
	# --remove_subgraphs 10 \
	# --output_modularity \
	# --output_graph_info \
	# --output_pickle \
	# --seq_to_ref /data/araju/ptpsinpdb/FinalProcessing/Final/3EAX_qFit_chainA.pdb \
	# -p

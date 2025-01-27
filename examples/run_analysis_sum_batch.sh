inputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run40_shivani_samples/'
outputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run40_shivani_samples/SumNetwork_Removed95_CommunitiesAll/'
outputfilename=${outputdir}'/SumNetwork'

mkdir -p ${outputdir}

pipenv run python ../scripts/analysis_sum.py \
	${inputdir}All_PTPsMultiNetwork.pkl \
	${outputfilename} \
	--remove_weak_edges 95 \
	--seq_to_ref /data/araju/ptpsinpdb/FinalProcessing/Final/1SUG_qFit_chainA.pdb \
	--remove_subgraphs 5 \
	--detect_communities \
	--output_modularity \
	--output_graph_info \
	 --output_pickle \
	# --n_communities 6
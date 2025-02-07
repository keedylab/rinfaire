inputdir='/Volumes/Akshay/ptpsinpdb/MultiRIN_Outputs/batchruns/run39/'
outputdir='/Volumes/Akshay/ptpsinpdb/MultiRIN_Outputs/batchruns/run39/SumNetwork_Removed95_CommunitiesAll/'
outputfilename=${outputdir}'/SumNetwork'

mkdir -p ${outputdir}

pipenv run python ../analysis_sum.py \
	${inputdir}MultiNetwork.pkl \
	${outputfilename} \
	--remove_weak_edges 95 \
	--seq_to_ref /Volumes/Akshay/ptpsinpdb/FinalProcessing/Final/1SUG_qFit_chainA.pdb \
	--remove_subgraphs 5 \
	--detect_communities \
	--output_modularity \
	--output_graph_info \
	 --output_pickle \
	# --n_communities 6
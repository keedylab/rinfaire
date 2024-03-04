
inputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run11/'
outputfilename='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run12/SumNetwork_Removed95'

pipenv run python ../scripts/analysis_sum.py \
	${inputdir}MultiNetwork.pkl \
	${outputfilename} \
	--remove_weak_edges 95 \
	--seq_to_ref /data/araju/ptpsinpdb/FinalProcessing/Final/6B8Z_qFit_chainA.pdb \
	-p \
	# --detect_communities \
	# --output_modularity \
	# --n_communities 6

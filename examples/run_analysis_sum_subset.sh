
inputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run13/'
outputfilename='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run13/SumNetwork_Subset_Removed50_STEP'

pipenv run python ../scripts/analysis_sum.py \
	${inputdir}MultiNetwork.pkl \
	${outputfilename} \
    --subset 'Protein_Name' 'STEP' \
	--remove_weak_edges 50 \
	--seq_to_ref /data/araju/ptpsinpdb/FinalProcessing/Final/6H8R_qFit_chainA.pdb \
	# --seq_to_ref /data/araju/ptpsinpdb/FinalProcessing/Final/6B8Z_qFit_chainA.pdb \
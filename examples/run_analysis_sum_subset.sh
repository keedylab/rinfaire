inputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run13/'
outputfilename='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run13/SumNetwork_Subset_Removed50_STEP'

pipenv run python ../scripts/analysis_sum.py \
	${inputdir}SubNetwork_All_PTP1B.pkl \
	${outputfilename} \
    --subset 'Ligand_Type' 'small_molecule' \
	--remove_weak_edges 50 \
	--seq_to_ref /data/araju/ptpsinpdb/FinalProcessing/Final/1SUG_qFit_chainA.pdb \
	# --seq_to_ref /data/araju/ptpsinpdb/FinalProcessing/Final/6B8Z_qFit_chainA.pdb \

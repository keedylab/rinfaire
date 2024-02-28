
inputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run11/'
# outputfilename='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run11/FoxComparison_Removed95'
# inputsetpath='/data/araju/ptpsinpdb/MultiRIN_Outputs/inputfiles/FoxSectors.csv'

# outputfilename='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run11/FoxClinMuts_Removed95'
# inputsetpath='/data/araju/ptpsinpdb/MultiRIN_Outputs/inputfiles/FoxClinMuts.csv'

outputfilename='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run11/FoxExperimentalMuts_Removed95'
inputsetpath='/data/araju/ptpsinpdb/MultiRIN_Outputs/inputfiles/FoxExperimentalMuts.csv'

pipenv run python ../scripts/analysis_residues_of_interest.py \
	${inputdir}SumNetwork_Removed95.pkl \
	${outputfilename} \
	--input_set ${inputsetpath} \
	--include_adjacent_residues /data/araju/ptpsinpdb/FinalProcessing/Final/6B8Z_qFit_chainA.pdb \
	--col All_Activity_Mutants
	#--col Clin_Muts
	#--col Sector_B


inputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run11/'
outputfilename='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run11/FoxComparison'
inputsetpath='/data/araju/ptpsinpdb/MultiRIN_Outputs/inputfiles/FoxSectors.csv'

pipenv run python ../scripts/analysis_residues_of_interest.py \
	${inputdir}SumNetwork_Removed95.pkl \
	${outputfilename} \
	--input_set ${inputsetpath} \



inputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run7/'
outputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run7/'

pipenv run python ../scripts/analysis_covariance.py \
	${inputdir}MultiNetwork.pkl \
	${outputdir} \
	--remove_weak_edges 0 \
	-c \
	-l DBSCAN

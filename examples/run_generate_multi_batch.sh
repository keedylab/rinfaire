
inputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/inputfiles/'
outputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run7/'

pipenv run python ../scripts/generate_multi.py \
	${inputdir}InputPDBs_V2.txt \
	${inputdir}PROMALS3D_PTPsAlignment_withoutArchael.fa \
	${outputdir}


inputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/inputfiles/'
outputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/test3/'

pipenv run python multirin/main_multi.py \
	${inputdir}InputPDBsTest.txt \
	${inputdir}PTPsAlignment_Promals3D.fa \
	${outputdir}


inputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/inputfiles/'
outputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/testruns/test18/'

pipenv run python ../scripts/generate_multi.py \
	${inputdir}InputPDBsTest.txt \
	${inputdir}PTPsAlignment_Promals3D.fa \
	${outputdir}


inputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/inputfiles/'
outputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/testruns/test21/'

pipenv run python ../scripts/generate_multi.py \
	${inputdir}InputPDBsTest.txt \
	${inputdir}PROMALS3D_PTPsAlignment_withoutArchael.fa \
	${outputdir} \
	--output_info \
	--add_adjacent_residues \

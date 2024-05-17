
inputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/inputfiles/'
outputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/testruns/test35/'

pipenv run python ../scripts/generate_multi.py \
	${inputdir}InputPDBsTest.txt \
	${inputdir}PROMALS3D_PTPsAlignment_withoutArchael.fa \
	${outputdir} \
	--metadata /data/araju/ptpsinpdb/MultiRIN_Outputs/inputfiles/PTPs_Metadata.csv \
	--output_info \
	--no_norm_struct

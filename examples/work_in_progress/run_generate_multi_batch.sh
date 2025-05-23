inputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/inputfiles/'
outputdir='/data/araju/ptpsinpdb/MultiRIN_Outputs/batchruns/run40/'

pipenv run python ../scripts/generate_multi.py \
	${inputdir}InputPDBs_V2.txt \
	${inputdir}PROMALS3D_PTPsAlignment_withoutArchael.fa \
	${outputdir} \
	--metadata /data/araju/ptpsinpdb/MultiRIN_Outputs/inputfiles/PTPs_Metadata.csv \
	--output_info \
	--norm_type log
	# --no_norm_struct

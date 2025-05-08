inputdir='/your/input/dir/'
outputfilename='/your/output/dir/output_filename’
outputdir='/your/output/dir/'

pipenv run python /multirin/installation/dir/generate_multi.py \
	${inputdir}All_PTPs_input.txt \
	${inputdir}PROMALS3D_PTPsAlignment.fa \
	${outputdir} \
	--metadata /multirin/installation/dir/examples/PTPs_Metadata.csv \
	--output_info \
	--norm_type log

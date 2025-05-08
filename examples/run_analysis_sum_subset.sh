inputdir='/your/input/dir/'
outputfilename='/your/output/dir/output_filename’

pipenv run python /multirin/installation/dir/analysis_sum.py \
	${inputdir}All_PTPs_MultiNetwork.pkl \
	${outputfilename} \
    	--subset 'Collection_Temperature' '100' \
	--remove_weak_edges 50 \
	--seq_to_ref reference.pdb \

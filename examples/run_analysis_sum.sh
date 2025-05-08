inputdir='/your/input/dir/'
outputfilename='/your/output/dir/output_filename'

pipenv run python /multirin/installation/dir/analysis_sum.py \
	${inputdir}All_PTPs_MultiNetwork.pkl \
	${outputfilename} \
	--remove_weak_edges 0 \
	--scale_sum_network struct \
	--output_graph_info \

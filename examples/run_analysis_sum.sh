input_dir="/your/output/directory/test_run_1/"
output_filename="${input_dir}/All_PTPs_SumNetwork"
output_dir="/your/output/directory/"
rinfaire_install_dir="/your/install/rinfaire/"
zenodo_dir="/your/download/zenodo_files/"

python ${rinfaire_install_dir}/analysis_sum.py \
        ${input_dir}/MultiNetwork.pkl \
        ${output_filename} \
        --seq_to_ref ${zenodo_dir}/qFit_Models_CatDom/1SUG_qFit_chainA.pdb \
        --output_graph_info \
        --output_pickle

# Optional
# --remove_weak_edges 95
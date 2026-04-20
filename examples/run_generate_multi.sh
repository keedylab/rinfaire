output_dir="/your/output/directory/test_run_1/"
rinfaire_install_dir="/your/install/rinfaire/"
zenodo_dir="/your/download/zenodo_files/"

mkdir -p ${output_dir}

python ${rinfaire_install_dir}/generate_multi.py \
        ${zenodo_dir}/All_PTPs_input.txt \
        ${zenodo_dir}/PROMALS3D_All_PTPs_aligned.fa \
        ${output_dir}
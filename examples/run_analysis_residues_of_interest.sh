input_dir="/your/output/directory/test_run_1/"
output_filename="${input_dir}/All_PTPs_SectorA"
rinfaire_install_dir="/your/install/rinfaire/"
zenodo_dir="/your/download/zenodo_files/"

python ${rinfaire_install_dir}/analysis_residues_of_interest.py \
        ${input_dir}/All_PTPs_SumNetwork_Rem95_All.pkl \
        ${output_filename} \
        --input_set ${zenodo_dir}/ResiduesOfInterest_Fig6.csv \
        --col 6b_Sector_A \
        --find_significance ${zenodo_dir}/qFit_Models_CatDom/1SUG_qFit_chainA.pdb \
        --n_iter_sig_test 100 \
        --cumulative_histogram \
        --no_normalize_by_total \
        --output_pickle

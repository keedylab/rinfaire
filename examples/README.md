# RINFAIRE Example Files & Scripts

## Input Files

Download files from our Zenodo:
```bash
curl --output zenodo_files.zip https://zenodo.org/api/records/18689236/files-archive 
unzip -d zenodo_files zenodo_files.zip
cd zenodo_files
unzip qFit_Models_FullStruct.zip
unzip qFit_Models_CatDom.zip
ls qFit_Models_CatDom/*.pdb > All_PTPs_input.txt
```

File Details:
- `All_PTPs_input.txt`: A `.txt` file containing the absolute paths to your qFit PDB files.
- `PROMALS3D_PTPsAlignment.fa`: A Multiple Sequence Alignment (MSA) file (`.fa`) of your structures
- `PTPs_Metadata.csv`: Optional metadata file in csv format to generate networks from subsets of structures

These files can be fed directly into the example scripts (just provide the path to your Zenodo download directory)

## Example Scripts

Below are example scripts for each command that can be run and compared with output files from our analyses present in the Zenodo and in the paper.

### 1. Generate a "multinetwork" containing individual networks for a set of structures using `generate_multi.py`:

```bash
bash run_generate_multi.sh
```

Outputs in Zenodo:
- All_PTPs_MultiNetwork.pkl

### 2. Generate a sum network for the entire multinetwork, including scaling and visualization using `analysis_sum.py`:

```bash
bash run_analysis_sum.sh
```

Note: To keep only top 5% of edges (for example scripts below) add the following flag (```--remove_weak_edges 95```)

Outputs in Zenodo:
- All_PTPs_SumNetwork.pkl
- All_PTPs_SumNetwork_Degree+Edges.xlsx

### 3. Generate a sum network for a **subset** of the multinetwork, including scaling and visualization *(alternative to step 2)* using `analysis_sum.py`:

```bash
bash run_analysis_sum_subset.sh
```

### 4. Evaluate overlap between sum network and user-defined residues of interest, including statistical significance:

```bash
bash run_analysis_residues_of_interest.sh
```

Outputs:
- Compare Jaccard index to those reported in Fig 6 & S16 (note that you will have to use top 5% sum network as input)

## Other

*Please note* that additional example scripts in the directory 'work in progress' may not have been tested with the current code so should be used with caution. 

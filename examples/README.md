# RINFAIRE

RINFAIRE (**R**esidue **I**nteraction **N**etworks **F**rom **A**lternate conformations **I**n **RE**lated structures) uses alternate conformations in crystallographic multiconformer models to calculate residue interaction networks (RINs) for sets of protein structures.  
It can then normalize and sum these networks, calculate differences between networks, and identify communities within a network, among other downstream analyses.

If you use this software, please cite:

[Raju A, Sharma S, Riley BT, Tan Y, Kim M, & Keedy DA. Mapping allosteric rewiring in related protein structures from collections of crystallographic multiconformer models. Preprint forthcoming (2025).](https://keedylab.org/publications/)

If you also use qFit, please cite:

[Wankowicz SA, Ravikumar A, Sharma S, Riley BT, Raju A, Hogan DW, van den Bedem H, Keedy DA, & Fraser JS. Automated multiconformer model building for X-ray crystallography and cryo-EM. eLife (2024).](https://elifesciences.org/articles/90606)

---

## Installation

`git clone -b main git@github.com:keedylab/RINFAIRE.git`

## Usage

RINFAIRE uses a series of commands to generate a sum network from a series of input structures and perform additional analyses. 
Example scripts for key steps/modules of RINFAIRE are listed in the [examples](examples/) directory.
The general steps of a typical workflow are as follows:

1. Preprocess the structures with [qFit](https://github.com/ExcitedStates/qfit-3.0) to accurately model any missing alternate conformations.
List the full paths of the qFit output models in a text file for input to RINFAIRE; see [examples](examples/).

2. If you wish to generate networks from subsets of structures, create a metadata file in csv format. 
See the metadata file in [examples](examples/); also supplied with the paper.

3. Perform a multiple sequence alignment for all input qFit models, using a tool such as [PROMALS3D](https://bio.tools/promals3d).
For more specific instructions, see the Methods section of the paper.

4. Generate a "multinetwork" consisting of aligned RINs for all structures:

`python scripts/generate_multi.py input_pdbs.txt /path/to/align.fa /output_dir/output_name --metadata /path/to/metadata.csv --output_info`

5. Normalize and sum the individual networks from the multinetwork to generate a sum network:

`python scripts/analysis_sum.py /input_dir/MultiNetwork.pkl /output_dir/output_name --output_graph_info`

6. Perform additional network analyses as desired; see the [examples](examples/) directory.

## License

The code is licensed under the MIT licence (see LICENSE).

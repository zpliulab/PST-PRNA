## Description
PST-PRNA is a method to decipher RNA binding sites on protein surface based on protein surface topography.
To achieve this, PST-PRNA builds the topographies and applies deep learning methods to learn from these.
For convenient use, please visit the web service www.zpliulab.cn/PSTPRNA.
For the standalone offline version, please install and use as follows.


## Standard alone Software prerequisites
* [Conda](https://docs.conda.io/en/latest/miniconda.html) Conda is recommended for environment management.
* [Python](https://www.python.org/) (3.6).
* [reduce](http://kinemage.biochem.duke.edu/software/reduce.php) (3.23). To add protons to proteins.
* [MSMS](http://mgltools.scripps.edu/packages/MSMS/) (2.6.1). To compute the surface of proteins.
* [PSI-Blast](https://blast.ncbi.nlm.nih.gov/) (2.6.0+) To generate PSSM.
* [hhblits](https://github.com/soedinglab/hh-suite) (3.3.0) To generate HMM.
* [CD-HIT](http://weizhong-lab.ucsd.edu/cd-hit/)(4.8.1) To clustering protein sequences .
* [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/) To standardize secondary structure assignment.

## Some important Python packages
* [BioPython](https://github.com/biopython/biopython) (1.78). To parse PDB files.
* [Pytorch](https://pytorch.org/) (1.7.1). pytorch with GPU version. Use to model, train, and evaluate the actual neural networks.
* [scikit-learn](https://scikit-learn.org/)(0.24.1).

## Specific usage
ONGOING

#1 Download and install the standard alone software listed above.
    Change the paths of these executable file at default_config/bin_path.py.


#2 Topography preparing
  a. The script 'protein.py' contains the class RBP which interates all procedures that are needed to convert a RBP to topographies.
  b. And for each protein, it takes tens of minutes to calculates topographies. So we recommend using parallel computing tools, such as [slurm](https://slurm.schedmd.com/).
    The bash script 'prepare_all.slurm' helps for extracting topographies in parallel cooperating with the python script 'prepare_all.py'.
  c. Users can also use 'prepare_all.py' all alone for preprocessing data. The files containing RBP_ids are in data/pdbid_chain. And the path of PDB_id lists should be specific
    the two 'prepare_all' scripts.

#2 Training
  To train an ab initio model, simply uses the script 'train.py'. Specific the RBPs list in default_config/dir_options:
  python train.py

#4 Predicting
  To predict new RNA-binding sites on newly protein, the topography preparing process is needed. Then
  uses the script 'predict.py':
  python predict.py

## Lisence
PST-PRNA is released under an [MIT Lisense](LICENSE)
## Reference
If you use this code, please use the bibtex entry in [citation.bib](citation.bib)
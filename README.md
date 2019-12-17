# SPARTA - SPARse peaks impuTAtion

SPARTA is a method for SPARse peaks impuTAtion that leverages predictive information within bulk ENCODE data to impute missing regions for a histone mark or transcription factor of interest. SPARTA is tested on a recent dataset (Grosselin et al. 2019) to impute DNA regions of H3K4me3 and H3K27me3 mark in B-cell and T-cell.

<img src="workflow/SPARTA.png" width="800">

## Installation with ANACONDA  

Install first anaconda in case you do not have it for linux machine. In any way we highle recommend to install the most recent one.

```
wget https://repo.anaconda.com/archive/Anaconda3-2019.10-Linux-x86_64.sh
bash Anaconda3-2019.10-Linux-x86_64.sh
```

Create a conda environment `sparta` with anaconda:

```
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda create -n sparta python=3.7 anaconda mpi4py
```
Finally activate the environment before running the algorithm.

`conda activate sparta`


## Clone repository and display command line arguments

```
git clone https://github.com/salbrec/SPARTA.git
cd SPARTA/

python SPART.py --help

usage: SPARTA.py [-h] --bed BED --targets TARGETS [--outdir OUTDIR]
                 [--genome GENOME] [--binsize BINSIZE]
                 [--estimators ESTIMATORS] [--simulate]

SPARTA: SPARse peaks impuTAtion

optional arguments:
  -h, --help            show this help message and exit
  --bed BED, -b BED     Sparse single-cell input dataset path
  --targets TARGETS, -t TARGETS
                        Target(s) defining the bulk reference data (ususally
                        the on used in scChIP). When multiple targets
                        provided, separate by "+"
  --outdir OUTDIR, -o OUTDIR
                        Output directory. Default: ./
  --genome GENOME, -g GENOME
                        Genome assembly
  --binsize BINSIZE, -bs BINSIZE
                        Size of the bins (genomic regions)
  --estimators ESTIMATORS, -e ESTIMATORS
                        Number of trees in Random Forest
  --simulate            Impute only 100 bins for testing the software


```

## Running SPARTA

Applying SPARTA on a bed file that contains the sparse input from an individual single cell is as simple as this:

```
python SPARTA.py -b ./scExamples/H3K4me3_hg38_5kb/BC9244609.bed -t H3K4me3
```
In order to test whether everyhing is running smoothly you can add `--simulate` that restricts SPARTA to train machine learning models for only 100 bins.

## SPARTA used with some arguments

Given an example from H3K27me3 in another resolution, the call could look like this, expecting the output to be stored in `/BC20160289/`.

`python SPARTA.py -b ./scExamples/H3K27me3_hg38_50kb/BC20160289.bed -t H3K27me3 -bs 50kb -o ./BC20160289/`

## Runtime and MPI

Due to the amount of machine learning models to be trained SPARTA can take up to 20h for the imputation of one cell for H3K4me3 in 5kb resolution. Though, SPARTA can be applied using MPI which we highly recommend. An example `slurm` script is provided that allows to run SPARTA using many cores from several compute nodes.  

## SPARTA output

Blaaa

#### Option 1 command line:

- imputation of a single cell from the command line 289 bins in 1 single cell BC9244609.bed:  

```
python sparta.py sc_examples/H3K4me3_hg38_5kb/BC9244609.bed H3K4me3 hg38 5kb RFC-100 -test
```

Check the `tes_out` folder where you will have the output of sparta corresponding to two files of the same size. One file in `.bed` format can be used for enrichment analysis, the second file in `.sparta` format has the information about:  

- the bin ID corresponding to the scChIP-seq region in the reference ENCODE portal    
- chr  
- start  
- end  
- frequency of the bin in the reference ENCODE dataset for the mark investigated (H3K4me3 in this case)
- imputed probability  

Please note that in `.sparta` format file, the bin already present in the single cell ChIP-seq file will have in the 6th column -1 (because imputation is not needed for that region). Furthermore in this simple test in the command line we will try to impute only 43 regions.    

#### Option 2 cluster: 

If you want to run sparta for all the regions in the ENCODE reference of approximately ~300000 regions we propose the `slurm` script that can run in a slurm scheduler. For this run:  

```
sbatch slurm
```
Within the slurm script the option `-test` is deactivated and all the bins present in the training test will be tested by SPARTA for classification.  

Sparta in this case will output the top ~30000 regions with the best imputed probability. We selected ~30000 regions because for the mark H3K4me3 is the average number in a bulk experiment.  

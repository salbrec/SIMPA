# SPARTAN - SPARse peaks impuTAtioN

SPARTAN is a method for SPARse peaks impuTAtioN that leverages predictive information within epigenomic data from ENCODE to impute missing protein-DNA interactions for a histone mark or transcription factor of interest. SPARTAN was tested on a recent dataset (Grosselin et al. 2019) to impute missing regions in sparse data from single-cell ChIP-seq of H3K4me3 and H3K27me3 in B-cells and T-cells. Different to common single-cell imputation methods, SPARTAN leverages predictive information within bulk ChIP-seq experimental data. This dataset contains > 2.200 experiments downloaded from ENCODE and available in this repository, preprocessed for SPARTAN. The user provides peaks of one single cell that are used by SPARTAN to impute missing interactions for the target (histone mark or transcription factor, also specified by the user) of interest while keeping cell-type specificity and the cells individuality. In our preprint on [bioRxiv](https://www.biorxiv.org/content/10.1101/2019.12.20.883983v1) we present SPARTAN's capability to complete sparse single-cell input while improving cell-type clustering and recovering cell-type-specific pathways.

<img src="figure/SPARTAN.png" width="900">

## Installation with ANACONDA  

First, install anaconda in case you do not have it in your linux machine. We highly recommend to install the most recent one.

```
wget https://repo.anaconda.com/archive/Anaconda3-2019.10-Linux-x86_64.sh
bash Anaconda3-2019.10-Linux-x86_64.sh

```
Accept licence and installation requirements with "return" and "yes", but follow the instructions, you might like to change the directory for anaconda. After installation it is necessary to initialize conda with:
```
source ~/.bashrc
```

### Create a conda environment `spartan` with anaconda:

```
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda create -n spartan python=3.7 anaconda pandas=0.25.1 numpy=1.17.2 mpi4py=3.0.2
```
Finally activate the environment before running the algorithm:

`conda activate spartan`


## Clone repository and display command line arguments

```
git clone https://github.com/salbrec/SPARTAN.git
cd SPARTAN/

python SPARTAN.py --help
```

The text of the help has to be as follow:

```
usage: SPARTAN.py [-h] --bed BED --targets TARGETS [--outdir OUTDIR]
                 [--genome {hg38}] [--binsize {5kb,50kb}]
                 [--estimators ESTIMATORS] [--simulate]

SPARTAN: SPARse peaks impuTAtioN

optional arguments:
  -h, --help            show this help message and exit
  --bed BED, -b BED     Sparse single-cell input dataset path
  --targets TARGETS, -t TARGETS
                        Target(s) defining the bulk reference data (ususally
                        the on used in scChIP). When multiple targets
                        provided, separate by "+"
  --outdir OUTDIR, -o OUTDIR
                        Output directory. Default: ./
  --genome {hg38}, -g {hg38}
                        Genome assembly
  --binsize {5kb,50kb}, -bs {5kb,50kb}
                        Size of the bins (genomic regions)
  --estimators ESTIMATORS, -e ESTIMATORS
                        Number of trees in Random Forest
  --simulate            Impute only 100 bins for testing the software

```

## Running SPARTAN

Applying SPARTAN on a bed file that contains the sparse input (300 peaks after removing gender specific chromosomes) from an individual single cell is as simple as this:

```
python SPARTAN.py -b ./scExamples/H3K4me3_hg38_5kb/BC8791969.bed -t H3K4me3 --simulate
```
The expected output looks like this:

```

Given the sparse input there are 300 genomic regions converted into 300 bins of size 5kb
Number of available bulk reference experiments: 178 (for H3K4me3)
Number of bins with a signal in at least one bulk: 314856
Number of candidate bins: 100 (simulation mode is on)
Shape of matrix for training features: (178, 300)

##### Pre-Processing is done ... #####

Training 92 models for 100 candidate bins
Random Forest is used with 100 trees
The task is shared by 1 processors

##### Writing output to "./" #####

Reference bulk experiments have in average 32584 bins
32284 bins were imputed (not really, because of "simulate")
Done!


```

Note, the optinal argument `--simulate` is used here for testnig purpose as it restricts SPARTAN to train machine learning models for only 100 candidate bins. To get the full result, remove it.

### Command line arguments

Given an example from H3K27me3 in another resolution, the call could look like this, expecting the output to be stored in `/BC20160289/`.  

```
python SPARTAN.py -b ./scExamples/H3K27me3_hg38_50kb/BC20160289.bed -t H3K27me3 -bs 50kb -o ./BC20160289/ --simulate
```

## SPARTAN output

SPARTAN creates two files reusing the filename of the given bed file as prefix. 

One in `spartan` format, a table describing one bin per line with the following columns:

- ID of the bin
- chromosome
- start 
- end
- reference frequency
- imputed probability

The bins on top of this table have no imputed probability as those are the bins observed for the given single cell. The following imputed bins are ranked by the imputed probability.

The second file is in `bed` format and describes the genomic regions for the single-cell bins followed by bins ranked by their imputed probability. The number of bins within this file is the average number of bins that SPARTAN calculates based on the target-specific reference experiments used within the algorithm.

It is up to the user to create further bed-files of different sizes derived from the SPARTAN table:

```
head -n <number of bins> IMPUTATION_freq_prob.spartan | awk '{print $2"\t"$3"\t"$4}' > new.bed
```

### Runtime and MPI

Due to the amount of machine learning models to be trained, SPARTAN can take up to 13h for the imputation of one cell in 5kb resolution for H3K4me3. However, SPARTAN is an MPI implementation that automatically distributes the computationally heavy part to multiple cores on a system providing an MPI installation. The number of cores can be defined by the user. An example `slurm` script is provided that allows to run SPARTAN using many cores from several compute nodes after adapting the account, partition, and other system-specific calls like "module load". We achieved the best efficiency by using a whole node (40 cores or an Intel® Xeon® Prozessor E5-2630 v4) that imputed one single cell within  minutes. However, the script provided uses 4 nodes with 160 cores doing the same job in ~8 minutes.

```
sbatch slurm
```
Having an MPI installation on a Linux server or local Linux machine, SPARTAN can also be used with `mpiexec`. The following call runs with 2 cores:
```
mpiexec -n 2 python SPARTAN.py -b ./scExamples/H3K4me3_hg38_5kb/BC8791969.bed -t H3K4me3
```

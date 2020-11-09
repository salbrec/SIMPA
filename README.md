# SIMPA - Single-cell chIp-seq iMPutAtion

SIMPA is a method for single-cell ChIP-seq imputation that leverages predictive information within epigenomic data from ENCODE to impute missing protein-DNA interactions for a histone mark or transcription factor of interest. SIMPA was tested on a recent dataset (Grosselin et al. 2019) to impute missing regions in sparse data from single-cell ChIP-seq of H3K4me3 and H3K27me3 in B-cells and T-cells. Different to common single-cell imputation methods, SIMPA leverages predictive information within bulk ChIP-seq experimental data. This dataset contains > 2.200 experiments downloaded from ENCODE and available in this repository, preprocessed for SIMPA. The user provides peaks of one single cell that are used by SIMPA to impute missing interactions for the target (histone mark or transcription factor, also specified by the user) of interest while keeping cell-type specificity and the cells individuality. In our preprint on [bioRxiv](https://www.biorxiv.org/content/10.1101/2019.12.20.883983v3) we present SIMPA's capability to complete sparse single-cell input while improving cell-type clustering and recovering cell-type-specific pathways.

<img src="figure/SIMPA.png" width="900">

## Installation

SIMPA, implemented in Python, runs on a Linux operating system and was tested on:
- Ubuntu 16.04.6 LTS (Xenial Xerus)
- Ubuntu 18.04.3 LTS (Bionic Beaver)
- CentOS Linux 7 (Core)
- Debian GNU/Linux 9 (stretch)
- macOS Catalina (10.15.2)

The following installation steps may take up to 15 minutes in total. Open a Linux terminal to run the commands for the installation and execution of SIMPA.

### Installation of ANACONDA  

First, install anaconda in case you do not have it in your linux machine. We highly recommend to install the most recent one.

```
wget https://repo.anaconda.com/archive/Anaconda3-2019.10-Linux-x86_64.sh
bash Anaconda3-2019.10-Linux-x86_64.sh

```
Accept licence and installation requirements with "return" and "yes", but follow the instructions, you might like to change the directory for anaconda. After installation it is necessary to initialize conda with:
```
source ~/.bashrc
```

### Create a conda environment `simpa` with anaconda:

```
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda create -n simpa python=3.7 anaconda pandas=0.25.1 numpy=1.17.2 mpi4py=3.0.2 terminaltables=3.1.0
```
Finally activate the environment before running the algorithm:

`conda activate simpa`


## Clone repository and display command line arguments
In case `git` is not installed, you can install it via `conda` as well:

```
conda install -c anaconda git
```

Clone repository and run SIMPA to get the usage information:

```
git clone https://github.com/salbrec/SIMPA.git
cd SIMPA/

python SIMPA.py --help
```

The expected output looks like this:

```
usage: SIMPA.py [-h] --bed BED --targets TARGETS [--outdir OUTDIR]
                [--genome {hg38}] [--binsize {5kb,50kb}]
                [--estimators ESTIMATORS] [--simulate]

SIMPA - Single-cell chIp-seq iMPutAtion

optional arguments:
  -h, --help            show this help message and exit
  --bed BED, -b BED     Path to bed file with sparse single-cell input
  --targets TARGETS, -t TARGETS
                        Target(s) defining the specific reference experiments
                        (ususally the one used in the scChIP). When multiple
                        targets are provided, separate by "+"
  --outdir OUTDIR, -o OUTDIR
                        Output directory. Default: "./"
  --genome {hg38}, -g {hg38}
                        Genome assembly
  --binsize {5kb,50kb}, -bs {5kb,50kb}
                        Size of the bins (genomic regions)
  --estimators ESTIMATORS, -e ESTIMATORS
                        Number of trees in Random Forest
  --simulate            Impute only 100 bins for testing the software

```

## Running SIMPA

Applying SIMPA on a bed file that contains the sparse input (300 peaks after removing gender specific chromosomes) from an individual single cell is as simple as this:

```
python SIMPA.py -b ./scExamples/H3K4me3_hg38_5kb/BC19017409_B-cell.bed -t H3K4me3 --simulate
```
The expected output looks like this:

```

Given the sparse input there are 643 genomic regions converted into 643 bins of size 5kb
Number of available bulk reference experiments: 178 (for H3K4me3)
Number of bins with a signal in at least one bulk: 314856
Number of candidate bins: 100 (simulation mode is on)
Shape of matrix for training features: (178, 643)

##### Pre-Processing is done ... #####

Training 92 models for 100 candidate bins
Random Forest is used with 100 trees
The task is shared by 1 processors

##### Writing output to "./" #####

Reference bulk experiments have in average 32584 bins
31941 bins were imputed (not really, because of "simulate")
Done!



```

Note, the optinal argument `--simulate` is used here for testnig purpose as it restricts SIMPA to train machine learning models for only 100 candidate bins. To get the full result, remove it.

### Command line arguments

Given an example from H3K27me3 in another resolution, the call could look like this, expecting the output to be stored in `/BC20160289/`.  

```
python SIMPA.py -b ./scExamples/H3K27me3_hg38_50kb/BC20160289.bed -t H3K27me3 -bs 50kb -o ./BC20160289/ --simulate
```

## SIMPA output

SIMPA creates two files reusing the filename of the given bed file as prefix. 

One in `simpa` format, a table describing one bin per line with the following columns:

- ID of the bin
- chromosome
- start 
- end
- reference frequency
- imputed probability

The bins on top of this table have no imputed probability as those are the bins observed for the given single cell. The following imputed bins are ranked by the imputed probability.

The second file is in `bed` format and describes the genomic regions for the single-cell bins followed by bins ranked by their imputed probability. The number of bins within this file is the average number of bins that SIMPA calculates based on the target-specific reference experiments used within the algorithm.

It is up to the user to create further bed-files of different sizes derived from the SIMPA table:

```
head -n <number of bins> IMPUTATION_freq_prob.simpa | awk '{print $2"\t"$3"\t"$4}' > new.bed
```

## InterSIMPA

The following line demonstrates an example call for InterSIMPA and the corresponding output. You'll see, InterSIMPA provides the results within seconds.

```
python InterSIMPA.py -b ./scExamples/H3K4me3_hg38_5kb/BC8791969_B-cell.bed -t H3K4me3 --summit chr19:35329500
```

This is the output:

```


Given the sparse input there are 300 genomic regions converted into 300 bins of size 5kb
Number of available bulk reference experiments: 178 (for H3K4me3)

Peak summit for given target region: chr19 35329500 35329501
Bin ID for given target region: 537938
Shape of matrix for training features: (178, 300)

##### Pre-Processing is done ... #####

Bin is present in given cell: NO
Reference Frequency of the given bin: 0.146
The imputed probability: 0.986


Features, represented by bins, most important for the model (importance > 1%)
additional annotations regarding the next gene (if bin overlaps gene-body, Distance==0)

+--------+------------+---------------------------+----------------+----------+----------+----------+----------------+-----------------------------------------------------+
| Bin ID | Importance | Genomic Region            | Next Gene (NG) | Dist TSS | Genebody | Dist TTS | NG Orientation | NG Description                                      |
+--------+------------+---------------------------+----------------+----------+----------+----------+----------------+-----------------------------------------------------+
| 122644 | 9.9%       | chr3:122075000-122080000  | CD86           | 22138    | overlap  | -43643   | plus           | CD86 molecule                                       |
| 219638 | 5.9%       | chr6:37005000-37010000    | FGD2           | 1845     | overlap  | -21572   | plus           | FYVE, RhoGEF and PH domain containing 2             |
| 523936 | 4.7%       | chr18:45685000-45690000   | LOC105372093   | 95336    | overlap  | -50450   | minus          | uncharacterized LOC105372093                        |
| 330217 | 4.4%       | chr9:114615000-114620000  | TMEM268        | 13603    | overlap  | -28922   | plus           | transmembrane protein 268                           |
| 574888 | 4.3%       | chr22:50315000-50320000   | DENND6B        | 9516     | overlap  | -8470    | minus          | DENN domain containing 6B                           |
| 245646 | 4.3%       | chr6:167045000-167050000  | FGFR1OP        | 48183    | overlap  | -5218    | plus           | FGFR1 oncogene partner                              |
| 168582 | 3.8%       | chr4:153470000-153475000  | TMEM131L       | 6154     | overlap  | -165108  | plus           | transmembrane 131 like                              |
| 94695  | 3.7%       | chr2:224520000-224525000  | CUL3           | 62897    | overlap  | -52350   | minus          | cullin 3                                            |
| 174767 | 3.5%       | chr4:184395000-184400000  | IRF2           | 77080    | overlap  | -14759   | minus          | interferon regulatory factor 2                      |
| 5263   | 3.4%       | chr1:26315000-26320000    | CD52           | -458     | overlap  | -3023    | plus           | CD52 molecule                                       |
| 563990 | 2.8%       | chr21:42530000-42535000   | SLC37A1        | 32878    | overlap  | -48936   | plus           | solute carrier family 37 member 1                   |
| 245477 | 1.9%       | chr6:166200000-166205000  | RNU6-153P      | -6903    | no       | -7006    | plus           | RNA, U6 small nuclear 153, pseudogene               |
| 214646 | 1.9%       | chr6:12045000-12050000    | HIVEP1         | 45404    | overlap  | -164534  | plus           | HIVEP zinc finger 1                                 |
| 324033 | 1.8%       | chr9:83695000-83700000    | UBQLN1         | 10479    | overlap  | -37532   | minus          | ubiquilin 1                                         |
| 137050 | 1.7%       | chr3:194105000-194110000  | LOC102724877   | 38184    | overlap  | -15974   | plus           | uncharacterized LOC102724877                        |
| 59334  | 1.6%       | chr2:47715000-47720000    | LOC100289315   | 12032    | no       | 11464    | plus           | NME/NM23 nucleoside diphosphate kinase 2 pseudogene |
| 124922 | 1.5%       | chr3:133465000-133470000  | BFSP2          | 69618    | overlap  | -7708    | plus           | beaded filament structural protein 2                |
| 257504 | 1.4%       | chr7:55530000-55535000    | VOPP1          | 40002    | overlap  | -98094   | minus          | VOPP1 WW domain binding protein                     |
| 390826 | 1.4%       | chr12:10390000-10395000   | KLRK1          | -2459    | overlap  | -20147   | minus          | killer cell lectin like receptor K1                 |
| 507990 | 1.3%       | chr17:49210000-49215000   | ABI3           | 2273     | overlap  | -10725   | plus           | ABI family member 3                                 |
| 399171 | 1.2%       | chr12:52115000-52120000   | LOC105369771   | 1539     | overlap  | -11789   | minus          | uncharacterized LOC105369771                        |
| 8077   | 1.2%       | chr1:40385000-40390000    | SMAP2          | 13828    | overlap  | -35826   | plus           | small ArfGAP2                                       |
| 387763 | 1.1%       | chr11:130160000-130165000 | ST14           | 2718     | overlap  | -47862   | plus           | suppression of tumorigenicity 14                    |
| 563856 | 1.1%       | chr21:41860000-41865000   | PRDM15         | 16982    | overlap  | -64275   | minus          | PR/SET domain 15                                    |
| 570912 | 1.0%       | chr22:30435000-30440000   | LOC105372990   | 4814     | overlap  | 999      | plus           | uncharacterized LOC105372990                        |
+--------+------------+---------------------------+----------------+----------+----------+----------+----------------+-----------------------------------------------------+



```

### Runtime and MPI

Due to the amount of machine learning models to be trained, SIMPA can take up to 5h for the imputation of one cell in 5kb resolution for H3K4me3. However, SIMPA is an MPI implementation that automatically distributes the computationally heavy part to multiple cores on a system providing an MPI installation. The number of cores can be defined by the user. An example `slurm` script is provided that allows to run SIMPA using many cores from several compute nodes after adapting the account, partition, and other system-specific calls like "module load". We achieved the best efficiency by using a whole node (40 cores or an Intel® Xeon® Prozessor E5-2630 v4) that imputed one single cell within 15 minutes. However, the script provided uses 4 nodes with 160 cores doing the same job in ~8 minutes.

```
sbatch slurm
```
Having an MPI installation on a Linux server or local Linux machine, SIMPA can also be used with `mpiexec`. The following call runs with 2 cores:
```
mpiexec -n 2 python SIMPA.py -b ./scExamples/H3K4me3_hg38_5kb/BC8791969.bed -t H3K4me3
```
Additional non-standard hardware for such an execution is not required.

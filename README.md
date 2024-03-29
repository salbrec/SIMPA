# SIMPA - Single-cell chIp-seq iMPutAtion

SIMPA is a method for **S**ingle-cell Ch**I**P-seq i**MP**ut**A**tion that leverages predictive information within epigenomic data from ENCODE to impute missing protein-DNA interactions for a histone mark or transcription factor of interest. SIMPA was tested on a recent dataset (Grosselin et al. 2019) to impute missing regions in sparse data from single-cell ChIP-seq of H3K4me3 and H3K27me3 in B-cells and T-cells. Different to common single-cell imputation methods, SIMPA leverages predictive information within bulk ChIP-seq experimental data. This dataset contains > 2.200 experiments downloaded from ENCODE and available in this repository, preprocessed for SIMPA. The user provides peaks of one single cell that are used by SIMPA to impute missing interactions for the target (histone mark or transcription factor, also specified by the user) of interest while keeping cell-type and single-cell specificity. In our paper published in [PLOS ONE](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0270043) we present SIMPA's capability to complete sparse single-cell input while maintaining cell-type clustering and recovering cell-type-specific pathways. In an additional analysis we demonstrate that imputed profiles are single-cell-specific and an application on single-cell data from mouse brain cells is included, as well. Moreover, SIMPA was extended by InterSIMPA which allows you to interpret the underlying machine learning models trained for the purpose of imputation. An example for how to do this, is provided below.

<img src="figure/SIMPA.png" width="900">

## Software Installation

SIMPA and InterSIMPA, both implemented in Python, run on a Linux operating system and was tested on:
- Ubuntu 16.04.6 LTS (Xenial Xerus)
- Ubuntu 18.04.3 LTS (Bionic Beaver)
- CentOS Linux 7 (Core)
- Debian GNU/Linux 9 (stretch)
- macOS Catalina (10.15.2)

The easiest way to get ready for SIMPA is pulling the docker and running the scripts inside the docker. The following descriptions explain how to get started with docker. However, it is also possible to install everything manually (see further installation guides at the bottom of this README). If you are familiar with the installation of python packages and you would simply like to update the installation you already have, these are the main packages needed together with the versions tested. Other versions might work perfectly fine, though.

```
pandas=1.1.3 
numpy=1.19.2 
mpi4py=3.0.3 
terminaltables=3.1.0
```

To start with docker, please open a Linux terminal and run the following commands to first install docker, then pull, and finally run the image.

```
sudo apt-get install docker
sudo apt-get install docker.io 

sudo docker login

sudo docker pull salbrec/simdock
sudo docker run -i -t -v "/home/:/home/" salbrec/simdock /bin/bash

```

When using the docker you’ll be directly within the right directory and able to use the scripts right away. Note, that you might save the output to your project directories and you might also prefer to use the updated github repository. Simply type `git pull` to update the repository within the docker.

We recommend to clone the git repository outside the docker. To clone the repository into your desired destination, change directory with `cd` and run:

```
git clone https://github.com/salbrec/SIMPA.git
```

#### *Check out further installation guides for running Docker with Windows 10 or creating your own conda environment (on the bottom of this README)*

Run SIMPA to get the usage information:

```
python SIMPA.py --help
```

The expected output looks like this:

```
usage: SIMPA.py [-h] --bed BED --targets TARGETS [--outdir OUTDIR] 
                [--genome {hg38,mm10}] [--binsize BINSIZE] 
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
  --genome {hg38,mm10}, -g {hg38,mm10}
                        Genome assembly
  --binsize BINSIZE, -bs BINSIZE
                        Size of the bins (genomic regions). For example "5kb" or "500bp"
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

The bins on top of this table have no imputed probability (value in imputed_probability column is set to -1) as those are the bins observed for the given single cell. The following imputed bins are ranked by the imputed probability in decreasing order (from 1.0 to 0.0).

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

## Additional parameters and options of InterSIMPA

By default, InterSIMPA lists all genes annotated by the single-cell regions with the highest feature importance values, but without any restriction for the distance to the gene. The parameter `–tssdist` can be used to exclude annotations with a distance to the TSS larger than the given value. Furthermore, InterSIMPA is able to include co-expression data from the STRING database as it was done for one of the validation analyses in our preprint. This functionality can be used by providing the Entrez symbol (gene name) of a gene that is related to the genomic position of interest (summit). The summit in the example above describes the center of the promoter for the gene CD22, one of the genes involved in the B-cell receptor signaling pathway. In order to extend the output of InterSIMPA the gene name has to be provided using the parameter `--gene` as in the following example:

```
python InterSIMPA.py -b ./scExamples/H3K4me3_hg38_5kb/BC8791969_B-cell.bed -t H3K4me3 --summit chr19:35329500 --gene CD22
```

The additional output will look like this:

```
You selected the gene "CD22" for an additional analysis based on the STRING co-expression data ...

+----------------------+--------------------+---+-----------------+---------------+
| InterSIMPA Relation  | Feature Importance |   | STRING Relation | Co-Expression |
+----------------------+--------------------+---+-----------------+---------------+
| CD86 -> CD22         | 9.93%              |   | CD22 <-> CD86   | 119           |
| FGD2 -> CD22         | 5.87%              |   | -               |               |
| LOC105372093 -> CD22 | 4.73%              |   | -               |               |
| TMEM268 -> CD22      | 4.42%              |   | -               |               |
| DENND6B -> CD22      | 4.33%              |   | -               |               |
| FGFR1OP -> CD22      | 4.26%              |   | -               |               |
| TMEM131L -> CD22     | 3.75%              |   | -               |               |
| CUL3 -> CD22         | 3.70%              |   | -               |               |
| IRF2 -> CD22         | 3.50%              |   | -               |               |
| CD52 -> CD22         | 3.42%              |   | CD22 <-> CD52   | 98            |
| SLC37A1 -> CD22      | 2.79%              |   | -               |               |
| RNU6-153P -> CD22    | 1.93%              |   | -               |               |
| HIVEP1 -> CD22       | 1.86%              |   | -               |               |
| UBQLN1 -> CD22       | 1.81%              |   | -               |               |
| LOC102724877 -> CD22 | 1.72%              |   | -               |               |
| LOC100289315 -> CD22 | 1.64%              |   | -               |               |
| BFSP2 -> CD22        | 1.47%              |   | -               |               |
| VOPP1 -> CD22        | 1.37%              |   | -               |               |
| KLRK1 -> CD22        | 1.36%              |   | CD22 <-> KLRK1  | 54            |
| ABI3 -> CD22         | 1.29%              |   | -               |               |
| LOC105369771 -> CD22 | 1.23%              |   | -               |               |
| SMAP2 -> CD22        | 1.16%              |   | -               |               |
| ST14 -> CD22         | 1.07%              |   | -               |               |
| PRDM15 -> CD22       | 1.06%              |   | -               |               |
| LOC105372990 -> CD22 | 1.04%              |   | -               |               |
+----------------------+--------------------+---+-----------------+---------------+
```

For most of the genes there is no co-expression data from STRING as indicated by the '-'. However, for the three genes that have a co-expression value one can clearly see that the higher the InterSIMPA feature importance, the stronger the co-expression.
Note that the corresponding STRING-based analysis in the preprint was done with hundreds of single-cell profiles which allows the collection of more genes for which both information are available, the InterSIMPA feature importance and the STRING co-expression. Accordingly, it was possible to compute correlation coefficients based on hundreds or even thousands of genes.

## Preparing bulk reference data for certain target(s) in a desired resolution

The script `prepareREFbins.py` implements the preprocessing for bulk reference experiments for any available target in a desired resolution. It can also be applied for several targets with one call when targets are separated by a '+'. The following example processes the reference experiments for H3K27ac and H3K9me3 in 2k resolution. Reference experiments are automatically downloaded as bed files and then preprocessed into bin sets. The resulting files are saved in a appropriate directory path and are then available for SIMPA or ItnerSIMPA.

```
python prepareREFbins.py -g hg38 -t H3K27ac+H3K9me3 -bs 2kb
```

### Runtime and MPI

Due to the amount of machine learning models to be trained, SIMPA can take up to 5h for the imputation of one cell in 5kb resolution for H3K4me3. However, SIMPA is an MPI implementation that automatically distributes the computationally heavy part to multiple cores on a system providing an MPI installation. The number of cores can be defined by the user. An example `slurm` script is provided that allows to run SIMPA using many cores from several compute nodes after adapting the account, partition, and other system-specific calls like "module load". We achieved the best efficiency by using a whole node (40 cores or an Intel® Xeon® Prozessor E5-2630 v4) that imputed one single cell within 15 minutes. However, the script provided uses 4 nodes with 160 cores doing the same job in ~8 minutes.

```
sbatch slurm
```

SIMPA can also be used with `mpiexec`. This works already within the docker or with a classical conda installation. The following call runs with, for instance, 2 cores:
```
mpiexec -n 2 python SIMPA.py -b ./scExamples/H3K4me3_hg38_5kb/BC8791969_B-cell.bed -t H3K4me3
```
Additional non-standard hardware for such an execution is not required.

## Further installation guides

### Installation with Docker Desktop on Windows

To install Docker Desktop follow the instructions on their website:
https://docs.docker.com/docker-for-windows/install/

Use git from powershell to clone seqQscorer
```
git clone https://github.com/salbrec/SIMPA.git

```
To get the image and activate it is similar to Linux. 
However it is advisable to only link the SIMPA folder.
Docker mentions that binding Windows Volumes can lead to performance drops and suggest to bind mounted folders from the linux filesystem in wsl rather than a folder on the windows file system.
Both work fine and can be accessed via powershell from the windows side or from the bash from the Linux/WSL side.

Below is an example from powershell, for linux just add sudo in front.
```
docker pull salbrec/simdock
docker run -i -t --name SIMPA -v "C:/Users/User/SIMPA:/SIMPA" salbrec/simdock 
```
Now you can just change to the newqscorer folder and start using the software!
```
(SIMPA) root@ xxx : cd SIMPA
```
In this example the SIMPA folder that is on windows is to find in the root of the docker image.
The docker image is named SIMPA and can be invoked by this name in the future.
You can copy files from the windows side and compute from the docker side.

Docker advises to use WSL, the mounted Linux System for Windows. If you want to use this, the installation and handling would be similar to the normal Linux installation, but the Installation of Docker Desktop for Windows also needs to be done.

### Installation with ANACONDA  

First, install anaconda in case you do not have it in your linux machine. We recommend to use the one that is suggested here. For the installation of Anaconda run the following two lines in your terminal.

```
wget https://repo.anaconda.com/archive/Anaconda3-2021.05-Linux-x86_64.sh
bash Anaconda3-2021.05-Linux-x86_64.sh

```
Accept licence and installation requirements with "return" and "yes", but follow the instructions, you might like to change the directory for anaconda. After installation it is necessary to initialize conda:
```
source ~/.bashrc

```

Now use the yml file `conda_env.yml` to create the conda environment and activate it.

```
conda env create -f conda_env.yml
conda activate simpa
```



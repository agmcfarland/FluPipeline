# FluPipeline

This pipeline takes a folder of illumina short read pairs in fastq.gz format and detects SNPs. 


## Installation instructions

Download the FluPipeline repository. The folder will be named `FluPipeline-main`.

Follow the directions below to create a conda environment with all required packages and dependencies.

Installing R packages will take a while.

```
# change directory to the pipeline folder
cd /path/to/FluPipeline-main

# create the conda environment
conda env create --name FluPipeline_env --file FluPipeline_env.yml

# enter the conda environment
conda activate FluPipeline_env


# enter R
R

# set the library that R will download packages to
rlib <- system('which R', intern=TRUE)
rlib <- sub(pattern='/bin/',replacement='/lib/',x=rlib)
rlib <- paste0(rlib,'/library')

#rlib = '/home/agmcfarland/miniconda3/envs/testenv/lib/R/library'

# verify that .libPaths() only includes the conda R library path
.libPaths(rlib)

# download R libraries. Do not update any libraries when prompted.
install.packages('remotes',repos='https://cloud.r-project.org/')
install_version("Rcpp", version = "1.0.7", repos = "http://cran.us.r-project.org", quiet=FALSE)
library(remotes)
install.packages('optparse',repos='https://cloud.r-project.org/', quiet=FALSE)
install.packages('optparse',repos='https://cloud.r-project.org/', quiet=FALSE)
install.packages('ggplot2',repos='https://cloud.r-project.org/', quiet=FALSE)
install.packages('knitr',repos='https://cloud.r-project.org/', quiet=FALSE)
install.packages('kableExtra',repos='https://cloud.r-project.org/', quiet=FALSE)
install.packages('stringr',repos='https://cloud.r-project.org/', quiet=FALSE)
install.packages('dplyr',repos='https://cloud.r-project.org/', quiet=FALSE)
install.packages('tidyverse',repos='https://cloud.r-project.org/', quiet=FALSE)
install_version('latticeExtra','0.6-28',repos='https://cloud.r-project.org/', quiet=FALSE) #2016 version
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite('ShortRead', suppressUpdates=TRUE, ask=FALSE)
biocLite('genbankr', suppressUpdates=TRUE, ask=FALSE)

q(save="no")

```


Afterwards change directory to the where the FluPipeline pipeline is located and type the following in the terminal:

```
python fluPipeline.py -h
```

If you got a help message, proceed to run the test below.


## Test

```
cd /path/to/FluPipeline-main

python fluPipeline.py \
--base_directory /path/to/output \
--sequence_directory /home/agmcfarland/flu_project/shared_data/test_data_6_samples \
--force \
--force_base_directory \
--threads 5

```

## Output files

Three main folders within the specified base_directory are created and contain the following.

1. sampleOutputs: Each sample read pair (from here on termed a sample) has a folder with the sample's name containing all files created during processing.

2. sampleResults: Each sample report is copied here for quick access.

3. sampleLogs: A record of how the run for each sample went.


Each run will produce a **run_summary.pdf**  report summarizing read coverage/depth, strain used as reference per sample, and whether errors occurred and for which samples an error occurred. 

Each run also will ouput a **runStats.csv** file.


## Usage

```
usage: fluPipeline.py [-h] [--base_directory BASE_DIRECTORY]
                      [--reference_directory REFERENCE_DIRECTORY]
                      [--sequence_directory SEQUENCE_DIRECTORY] [--force]
                      [--force_base_directory] [--threads THREADS]

optional arguments:
  -h, --help            show this help message and exit
  --base_directory BASE_DIRECTORY
                        directory that run samples will be saved in
  --reference_directory REFERENCE_DIRECTORY
                        directory containing reference strain files (.gb
                        format). Default is the references folder where fluPipeline.py is located
  --sequence_directory SEQUENCE_DIRECTORY
                        directory containing fastq sequence files (.gz format)
  --force               overwrite existing files
  --force_base_directory
                        overwrite existing directory
  --threads THREADS     number of processors to use for multiprocessing
  ```



# FOR DEVELOPMENT

```
# create a new conda environment call FluPipeLine_env
conda create --name FluPipeline_env3

# enter the environment
source activate FluPipeline_env #or conda activate FluPipeline_env


# download the following. press y when prompted
conda install -c conda-forge r=3.4.1

#'y'

conda install -c anaconda python=3.6.3

#'y'

pip install biopython numpy pandas

conda install -c bioconda bwa=0.7.15

#'y'

conda install -c bioconda samtools=1.7

#'y'

conda install -c bioconda bcftools=1.8

#'y'

conda install -c bioconda fastp=0.12.4

#'y'

conda install -c bioconda bbmap=38.18

#'y'

conda install -c anaconda ipython

#'y'

## enter R 

R

# download the following. Press an empty space when prompted.

install.packages('remotes',repos='https://cloud.r-project.org/', quiet=FALSE)
library(remotes)
install.packages('optparse',repos='https://cloud.r-project.org/', quiet=FALSE)
install.packages('ggplot2',repos='https://cloud.r-project.org/', quiet=FALSE)
install.packages('knitr',repos='https://cloud.r-project.org/', quiet=FALSE)
install.packages('kableExtra',repos='https://cloud.r-project.org/', quiet=FALSE)
install.packages('stringr',repos='https://cloud.r-project.org/', quiet=FALSE)
install.packages('dplyr',repos='https://cloud.r-project.org/', quiet=FALSE)
install.packages('tidyverse',repos='https://cloud.r-project.org/', quiet=FALSE)
install_version('latticeExtra','0.6-28',repos='https://cloud.r-project.org/', quiet=FALSE) #2016 version
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite('ShortRead', suppressUpdates=TRUE, ask=FALSE)
biocLite('genbankr', suppressUpdates=TRUE, ask=FALSE)

#press space

q(save="no")

```
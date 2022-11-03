#!/bin/bash

conda deactivate

conda remove -n testenv --all

yes | conda create -n testenv

conda activate testenv

conda config --append channels conda-forge
conda config --append channels bioconda
conda config --append channels main
conda config --append channels anaconda

yes | mamba install -c conda-forge gsl #needed to solve issuse with bcftools

yes | mamba install -c bioconda bcftools=1.15.1

yes | pip install biopython numpy pandas psutil pytz

yes | mamba install -c bioconda bwa=0.7.17

yes | mamba install -c bioconda bamutil=1.0.15

yes | mamba install -c bioconda fastp=0.22.0

yes | mamba install -c bioconda nextclade=2.4.0

yes | mamba install -c bioconda bbmap=38.18

# yes | mamba install -c bioconda spades=3.13.0

yes | mamba install -c bioconda lofreq=2.1.5

yes | mamba install -c bioconda samtools=1.15.1

yes | mamba install -c anaconda ipython

yes | mamba install -c bioconda freebayes=1.3.6

yes | mamba install -c bioconda picard-slim=2.27.4

# check install is successful

bcftools

samtools

lofreq

callvariants.sh

nextclade

bwa

freebayes

picard
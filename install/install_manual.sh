#!/bin/bash

yes | conda create --name testenv3

yes | conda activate testenv3

yes | conda install -c conda-forge r=3.4.1

yes | conda install -c bioconda bcftools=1.15 #maybe not needed

yes | conda install -c anaconda python=3.6.3

yes | pip install biopython numpy pandas

yes | conda install -c bioconda bwa=0.7.15

yes | conda install -c bioconda samtools=1.7

yes | conda install -c bioconda bamutil=1.0.15

yes | conda install -c bioconda fastp=0.22.0

yes | conda install -c bioconda bbmap=38.18

yes | conda install -c anaconda ipython

yes | Rscript ./r_packages.R
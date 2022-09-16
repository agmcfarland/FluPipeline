#!/bin/bash

conda deactivate

conda remove -n testenv --all

yes | conda create -n testenv

conda activate testenv

yes | pip install biopython numpy pandas psutil

yes | mamba install -c bioconda bwa=0.7.17

yes | mamba install -c bioconda bamutil=1.0.15

yes | mamba install -c bioconda fastp=0.22.0

yes | mamba install -c bioconda nextclade=2.4.0

yes | mamba install -c bioconda bbmap=38.18

yes | mamba install -c bioconda spades=3.13.0

yes | mamba install -c conda-forge libgcc-ng

yes | mamba install -c bioconda lofreq=2.1.5

yes | mamba install -c bioconda samtools=1.15.1

yes | mamba install -c anaconda ipython

# yes | mamba install -c bioconda lofreq


conda install -c conda-forge libzlib

conda install -c bioconda libdeflate=1.2

conda install -c conda-forge libgcc-ng=12.1.0

yes | mamba install -c bioconda bcftools=1.15.1




yes | mamba install -c bioconda samtools=1.15.1

yes | pip install biopython numpy pandas psutil

yes | mamba install -c bioconda bwa=0.7.17

yes | mamba install -c bioconda bamutil=1.0.15

yes | mamba install -c bioconda fastp=0.22.0

yes | mamba install -c bioconda nextclade=2.4.0

yes | mamba install -c bioconda bbmap=38.18


yes | mamba install -c anaconda ipython

# yes | mamba install -c bioconda lofreq

yes | mamba install -c bioconda bcftools=1.15.1

yes | mamba install -c bioconda spades=3.13.0

yes | mamba install -c conda-forge libgcc-ng

yes | mamba install -c bioconda lofreq=2.1.5

yes | mamba install -c bioconda samtools=1.15.1

yes | mamba install -c anaconda ipython

# yes | mamba install -c bioconda lofreq

yes | mamba install -c bioconda bcftools=1.15.1
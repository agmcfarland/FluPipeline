## instructions: when asked if you want to update any packages always say no (press `n`)

# installation
conda create --name FluPipeline_env

source activate FluPipeline_env

conda install -c conda-forge r=3.4.1

'y'

conda install -c anaconda python=3.6.3

'y'

pip install biopython numpy pandas

conda install -c bioconda bwa=0.7.15

'y'

conda install -c bioconda samtools=1.7

'y'

conda install -c bioconda bcftools=1.8

conda install -c bioconda fastp=0.12.4

'y'

conda install -c bioconda bbmap=38.18

'y'

conda install -c anaconda ipython

'y'

R

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

#0 or empty line?

q(save="no")
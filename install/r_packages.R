# set the library that R will download packages to
rlib <- system('which R', intern=TRUE)
rlib <- sub(pattern='/bin/',replacement='/lib/',x=rlib)
rlib <- paste0(rlib,'/library')
.libPaths(rlib)

# verify that .libPaths() only includes the conda R library path
.libPaths()

# download R libraries. Do not update any libraries when prompted.
install.packages('remotes',repos='https://cloud.r-project.org/')
library(remotes)
install_version('rlang', version = "0.4.12", repos = "http://cran.us.r-project.org", quiet=FALSE)
install_version('openssl', version = "1.4.5", repos = "http://cran.us.r-project.org", quiet=FALSE)
install_version("Rcpp", version = "1.0.7", repos = "http://cran.us.r-project.org", quiet=FALSE)
install.packages('optparse',repos='https://cloud.r-project.org/', quiet=FALSE)
install.packages('ggplot2',repos='https://cloud.r-project.org/', quiet=FALSE)
install.packages('knitr',repos='https://cloud.r-project.org/', quiet=FALSE)
install.packages('stringr',repos='https://cloud.r-project.org/', quiet=FALSE)
install.packages('dplyr',repos='https://cloud.r-project.org/', quiet=FALSE)
install.packages('tidyr',repos='https://cloud.r-project.org/', quiet=FALSE)
install_version('latticeExtra','0.6-28',repos='https://cloud.r-project.org/', quiet=FALSE) #2016 version
install.packages('lubridate',repos='https://cloud.r-project.org/', quiet=FALSE)
# install.packages('kableExtra',repos='https://cloud.r-project.org/', quiet=FALSE) # already gets installed?
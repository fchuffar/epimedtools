conda create --name epimedtools_env
conda activate epimedtools_env


conda install -n base conda-forge::mamba
mamba install -c anaconda -c bioconda -c conda-forge -c r jq
mamba install -c anaconda -c bioconda -c conda-forge -c r bioconductor-genomeinfodbdata 
mamba install -c anaconda -c bioconda -c conda-forge -c r r-base libopenblas bioconductor-geoquery bioconductor-affy bioconductor-biobase r-seqinr r-rcpparmadillo r-devtools r-fastmap r-matrix r-kernsmooth r-catools r-gtools r-nortest r-survival r-beanplot r-gplots bioconductor-hgu133plus2cdf

# under R
# devtools::install_github("fchuffar/epimedtools")

conda create --name epimedtools_env
conda activate epimedtools_env

conda install  -c anaconda libopenblas
conda install -c conda-forge r-base=4.2.0 r-devtools
conda install -c r r-rcpparmadillo 
conda install -c bioconda bioconductor-geoquery bioconductor-affy  bioconductor-biobase r-seqinr

# under R
# devtools::install_github("fchuffar/epimedtools")

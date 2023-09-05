conda create --name epimedtools_env
conda activate epimedtools_env

conda install  -c anaconda libopenblas
conda install -c conda-forge r-base=4.2.0
conda install -c bioconda bioconductor-geoquery bioconductor-affy  bioconductor-biobase r-seqinr
conda install -c r r-rcpparmadillo 
conda install -c conda-forge r-devtools

# under R
# devtools::install_github("fchuffar/epimedtools")

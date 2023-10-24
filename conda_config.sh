conda create --name epimedtools_env
conda activate epimedtools_env


conda install -n base conda-forge::mamba
mamba install -c anaconda -c bioconda -c conda-forge -c r libopenblas=0.3.21 r-base=4.2.0 bioconductor-geoquery=2.66.0 bioconductor-affy=1.76.0 bioconductor-biobase=2.58.0 r-seqinr=4.2_30 r-rcpparmadillo=0.12.6.4.0 r-devtools=2.4.5 r-fastmap=1.1.1 libstdcxx-devel_linux-64=11.2.0 r-matrix r-kernsmooth r-catools r-gtools r-nortest r-survival r-beanplot r-gplots

# under R
# devtools::install_github("fchuffar/epimedtools")

conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda
if test -e $HOME/miniconda3/envs/exfi_validation ; then
else
    conda env create --name exfi_validation --file requirements.txt
fi

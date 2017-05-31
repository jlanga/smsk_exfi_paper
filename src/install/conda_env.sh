conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda

if [ ! -f /tmp/foo.txt ]; then
    echo "File not found!"
fi
if [ ! -e $HOME/miniconda3/envs/exfi_validation ] ; then
    conda env create --name exfi_validation --file requirements.txt
else
    source activate exfi_validation
    conda install --yes --file requirements.txt
fi

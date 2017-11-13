if [ ! -e $(which conda)/.miniconda3/envs/exfi_validation ] ; then
    conda env create --name exfi_validation --file environment.yml
else
    source activate exfi_validation
    conda env update --file environment.yml
fi

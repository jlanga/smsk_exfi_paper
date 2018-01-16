#!/usr/bin/env bash

# Create conda environment
conda env update --file environment.yml
conda clean --all --yes

export PATH="$HOME/miniconda3_$TRAVIS_OS_NAME/bin:$PATH"
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda

source activate exfi_validation

pushd /opt/

pushd /opt/

# SDSL-lite
# https://hub.docker.com/r/adamnovak/sequence-graphs/~/dockerfile/
if [[ ! -d sdsl-lite/ ]]; then
    git clone https://github.com/simongog/sdsl-lite.git
fi
pushd sdsl-lite/ && \
sudo ./install.sh /usr/ && \
popd

# biobloomtools
# as in https://github.com/bcgsc/biobloom
if [[ ! -d biobloom/ ]]; then
    git clone https://github.com/bcgsc/biobloom.git
fi
pushd biobloom/ && \
git submodule update --init && \
./autogen.sh && \
./configure --prefix=/usr/local/ && \
make -j 2 && \
sudo make install && \
popd


# exfi package via github
git clone https://github.com/jlanga/exfi.git
pushd exfi/
python3 setup.py test
pip install . --user --no-deps --upgrade
popd

popd

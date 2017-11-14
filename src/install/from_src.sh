#!/usr/bin/env bash

mkdir -p src/
mkdir -p bin/

pushd src/

# SDSL-lite
# https://hub.docker.com/r/adamnovak/sequence-graphs/~/dockerfile/
if [[ ! -d sdsl-lite/ ]]; then
    git clone https://github.com/simongog/sdsl-lite.git
fi
pushd sdsl-lite/ && \
sudo ./install.sh /usr/local/ && \
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

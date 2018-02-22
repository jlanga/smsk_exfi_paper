#!/usr/bin/env bash

mkdir -p src/
mkdir -p bin/

pushd src/

# SDSL-lite
# https://hub.docker.com/r/adamnovak/sequence-graphs/~/dockerfile/
pushd sdsl-lite/ && \
sudo ./install.sh /usr/local/ && \
popd

# biobloomtools
# as in https://github.com/bcgsc/biobloom
pushd biobloom/ && \
git submodule update --init && \
./autogen.sh && \
./configure --prefix=/usr/local/ && \
make -j 2 && \
sudo make install && \
popd

popd

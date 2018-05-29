#!/usr/bin/env bash
set -euxo pipefail

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
git checkout 0a42916922d42611a087d4df871e424a8907896e && \
./autogen.sh && \
./configure --prefix=/usr/local/ && \
make -j 2 && \
sudo make install && \
popd

popd

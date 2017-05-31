#!/usr/bin/env bash

mkdir -p src/
mkdir -p bin/

pushd src/

# exfi package via github
git clone https://github.com/jlanga/exfi.git
pushd exfi/
python setup.py test
pip install . --no-deps --upgrade
popd



# biobloomtools via github
git clone https://github.com/bcgsc/biobloom biobloom
pushd biobloom/
./autogen.sh
./configure
make -j
cp BioBloomMaker/biobloommaker BioBloomCategorizer/biobloomcategorizer ../../bin/
popd

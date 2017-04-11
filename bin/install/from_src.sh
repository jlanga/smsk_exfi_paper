#!/usr/bin/env bash

mkdir -p src/

pushd src/

# exfi package via github
git clone https://github.com/jlanga/exon_finder.git .
pushd exon_finder/
python setup.py test
pip3 install --user --editable .
popd

popd
#!/usr/bin/env bash

mkdir -p src/

pushd src/

# exfi package via github
git clone https://github.com/jlanga/exfi.git
pushd exfi/
python setup.py test
pip install --editable .
popd

popd
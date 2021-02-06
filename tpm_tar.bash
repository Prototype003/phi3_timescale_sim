#!/bin/bash

# $1 is tpm type

pushd tpms/$1
find -type f -name '*.mat' > ls_out
tar -cvf -T ls_out $1.tar

# remove files
mv $1.tar ../ # Move tar out of directory
rm -f * # Empty the directory
mv ../$1.tar ./ # Move the tar back

# Pre-extract params.mat
tar -xf $1.tar params.mat

popd
#!/bin/bash

# DELETES all .mat files from split/$1
# $1 is tpm type

# Move joined result .mat file out
mv split/$1/joined.mat ./

# Remove .mat files from directory
pushd split/$1
rm -f *.mat
popd

# Return joined result .mat file
mv joined.mat split/$1/
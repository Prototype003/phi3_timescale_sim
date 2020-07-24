#!/bin/bash

# DELETES all .mat files from $1/
# $1 is tpm type

# Move params and joined .mat files out
mv $1/params.mat $1/tpms.mat ./

# Remove .mat files from directory
pushd $1
rm -f *.mat
popd

# Return joined result .mat file
mv params.mat tpms.mat $1/
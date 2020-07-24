#!/bin/bash

# Extracts all files from phis.tar in split/$1
# $1 is tpm type

pushd split/$1
tar -xf phis.tar
popd
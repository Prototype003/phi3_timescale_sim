#!/bin/bash

# Extracts all files from $1/$1.tar
# $1 is tpm type

pushd $1/
tar -xf $1.tar
popd
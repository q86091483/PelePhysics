#!/bin/bash

# Find first executable name available
make COMP=llvm TPL
make COMP=llvm
execname=`find . -name "Pele*.ex" | head -1`

if [ -z "$execname" ]
then
    make COMP=llvm TPL
    make COMP=llvm
    execname=`find . -name "Pele*.ex" | head -1`
fi


if [[ -f "PPreaction.txt" ]]; then
    rm PPreaction.txt
fi
if [[ -f "log" ]]; then
    rm log
fi
$execname input



python compareResults.py

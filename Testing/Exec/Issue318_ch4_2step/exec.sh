#!/bin/bash

# Compile
make COMP=llvm TPL
make COMP=llvm
execname=`find . -name "Pele*.ex" | head -1`

# Execute HR test
if [[ -f "PPreaction.txt" ]]; then
    rm PPreaction.txt
fi
$execname input



# Plot results
python compareResults.py

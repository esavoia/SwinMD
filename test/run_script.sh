#!/bin/bash

TARGET_EXE="SwinMD.exe"

echo
echo .. running SwinMD ...
echo

mpirun -n 4 ${TARGET_EXE} > stdout

# move outputs
mkdir -p outputs
mv *.out outputs/

echo "... Done! :D"

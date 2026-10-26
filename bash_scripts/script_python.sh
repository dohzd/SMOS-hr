#!/bin/bash

#OAR -l /nodes=1/core=16,walltime=48:0:0
#OAR --stdout tmp/script_python.%jobid%.out
#OAR --stderr tmp/script_python.%jobid%.err
#OAR --project snowem

ROOT_DIR=$1
DL_DIR=$2
AREA=$3
GRID=$4

export ROOT_DIR
echo "$ROOT_DIR"

export DL_DIR
echo "$DL_DIR"

export AREA
echo "$AREA"

export GRID
echo "$GRID"

export PATH="$HOME/.pixi/bin:$PATH"

cd "${ROOT_DIR}/prog" || exit 1

pixi run -e python python SMOS-hr-main/src/SMOS-hr/python/SMOSHR.py

echo "Finished"

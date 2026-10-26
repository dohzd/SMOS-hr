#!/bin/bash

#OAR -l /core=4,walltime=36:00:0
#OAR --stdout tmp/script_octave.%jobid%.out
#OAR --stderr tmp/script_octave.%jobid%.err
#OAR --project snowem
#####OAR -t devel

ROOT_DIR=$1
DL_DIR=$2
AREA=$3

if [[ "$AREA" == "Greenland" ]]; then
    reader="reader_Greenland.m"
fi

if [[ "$AREA" == "Antartica" ]]; then
    reader="reader_Antartica.m"
fi

export PATH="$HOME/.pixi/bin:$PATH"

# Pixi workspace located in this directory
cd "${ROOT_DIR}/prog"

# Search below a day, month or year directory
find "$DL_DIR" -type d -name 'SM*_1' -prune \
    -exec pixi run -e octave octave \
        --path SMOS-hr-main/src/SMOS-hr/octave \
        SMOS-hr-main/src/SMOS-hr/octave/${reader} {} +



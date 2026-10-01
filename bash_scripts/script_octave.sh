#!/bin/bash

#OAR -l /core=1,walltime=0:30:0
#OAR --stdout tmp/script_octave.%jobid%.out
#OAR --stderr tmp/script_octave.%jobid%.err
#OAR --project snowem
#####OAR -t devel

DIR=$1
# "/bettik/PROJECTS/pr-snowem/zeigerp/SMOS_data/L1C/2020/01"
source /applis/site/guix-start.sh
cd /home/zeigerp/programs/L1c_reader

PARAM=($DIR/SM*_1)

echo "$PARAM"
octave reader_Antartica.m $PARAM




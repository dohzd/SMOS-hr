#!/bin/bash

#OAR -l /nodes=1/core=32,walltime=0:30:00
#OAR --stdout tmp/script_octave.%jobid%.out
#OAR --stderr tmp/script_octave.%jobid%.err
#OAR --project snowem
#####OAR -t devel

DIR=$1
echo "$DIR"
DIR= "/bettik/PROJECTS/pr-snowem/zeigerp/SMOS_data/L1C/2020/01"
echo "$DIR"
source /applis/site/guix-start.sh
cd ~/programs/L1c_reader

#PARAM=($DIR/*/*/SM*_1)
#PARAM=($DIR/*/SM*_1)
PARAM=($DIR/*)


for PAR in ${PARAM[@]}
do
	cd /home/zeigerp/programs
	oarsub -S "./script_octave.sh $PAR"
done


# Old
#num_procs=32
#num_jobs="\j"  # The prompt escape for number of jobs currently running
#for PAR in ${PARAM[@]}
#do
#        while (( ${num_jobs@P} >= num_procs ))
#        do
#                wait -n
#        done
#        #octave reader_Greenland.m $PAR &
#        octave reader_Antartica.m $PAR &
#done


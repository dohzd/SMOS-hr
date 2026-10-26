#!/bin/bash

#OAR -l /nodes=1/core=2,walltime=2:00:00
#OAR --stdout tmp/script_octave.%jobid%.out
#OAR --stderr tmp/script_octave.%jobid%.err
#OAR --project snowem


DIR=$1
echo "$DIR"
DIR= "/silenus/PROJECTS/pr-snowem/huzardd/SMOS/data"
echo "$DIR"
source /applis/site/guix-start.sh
cd ~/programs/L1c_reader

#PARAM=($DIR/*/*/SM*_1)
#PARAM=($DIR/*/SM*_1)
PARAM=($DIR/*)


for PAR in ${PARAM[@]}
do
	cd /silenus/PROJECTS/pr-snowem/huzardd/prog/SMOS-hr-main/bash_scripts/
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


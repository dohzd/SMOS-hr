#!/bin/bash

#OAR -l /nodes=1/core=32,walltime=3:0:0
#OAR --stdout tmp/script_python.%jobid%.out
#OAR --stderr tmp/script_python.%jobid%.err
#OAR --project snowem

DIR=$1 #"/bettik/PROJECTS/pr-snowem/zeigerp/SMOS_data/L1C/2020/01"   #$1

source /applis/environments/conda.sh
#eval "$(conda shell.bash hook)"
#eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate smoshr

cd /home/zeigerp/programs/python
PARAM="$DIR"
export PARAM
echo "$PARAM"
NAME="Greenland"
export NAME
echo "$NAME"
python3 /home/zeigerp/programs/python/SMOSHR.py

#for PARAM in $DIR/*/
#do
#       echo "$PARAM"
#       export PARAM
#       python3 SMOSHR.py && sleep 5 & 
#       #/home/zeigerp/.conda/envs/smoshr/bin/conda run -n smoshr -p /home/zeigerp/.conda/envs/ /home/zeigerp/.conda/envs/smoshr/bin/python /home/zeigerp/programs/python/SMOSHR.py &
#       #/home/zeigerp/.conda/envs/smoshr/bin/python3 /home/zeigerp/programs/python/SMOSHR.py &
#       /home/zeigerp/.conda/envs/smoshr/bin/python "/home/zeigerp/programs/python/SMOSHR.py" &
#       /home/zeigerp/.conda/envs/smoshr/bin/python "print('Hello World')" &
#       #conda run -p "/home/zeigerp/.conda/envs/smoshr" python "/home/zeigerp/programs/python/SMOSHR.py" &
#done

echo "Finished"

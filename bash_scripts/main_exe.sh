#!/bin/bash

#OAR -l /core=1,walltime=12:00:00
#OAR --stdout tmp/main_exe.%jobid%.out
#OAR --stderr tmp/main_exe.%jobid%.err
#OAR --project snowem

while getopts "y:m:d:" param; do
    case "$param" in
        y) year="$OPTARG" ;;
        m) month="$OPTARG" ;;
        d) day="$OPTARG" ;;

    esac
done

if [[ -z "$year" || ( -z "$month" && -n "$day" ) ]]; then
    echo "Wrong parameter"
    exit 1
fi

date_path="$year"
if [[ -n "$month" ]]; then
    date_path="$date_path/$month"
fi
if [[ -n "$day" ]]; then
    date_path="$date_path/$day"
fi

name_area="Greenland"

grid="12"

url="ftps://smos-diss.eo.esa.int:990/SMOS/L1CL/MIR_SCLF1C/${date_path}/"

root_dir="/silenus/PROJECTS/pr-snowem/huzardd/SMOS/"
data_dir="${root_dir}data/L1C/"
download_dir="${data_dir}${date_path}"

# Wget reads the ESA login/password from ~/.netrc
wget -m -c --no-host-directories --cut-dirs=3 --ftps-implicit --no-check-certificate "$url" -P "$data_dir" || exit "$?"

find "$download_dir" -iname '*.zip' -type f -exec sh -c 'unzip -o -d "${0%.*}" "$0" && rm "$0"' '{}' ';'

# Capture the Octave job ID to make Python wait for this job
octave_submission=$(oarsub -S "./script_octave.sh $root_dir $download_dir $name_area") || exit 1
octave_job_id=$(echo "$octave_submission" | sed -n 's/^OAR_JOB_ID=//p')
# [[ -n "$octave_job_id" ]] || exit 1

oarsub -a "$octave_job_id" -S "./script_python.sh $root_dir $download_dir $name_area $grid"

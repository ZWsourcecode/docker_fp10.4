#!/bin/bash
### fix the missing footprint simulation and/or 14C concertration calculation
### script should be in parent folder of options
### run : bash fix_global.sh 11 inputs 
### zhendong.wu@nateko.lu.se

F_PATH="/flexpart"
FLEXPARTPATH="/usr/local/flexpart_v10.4/src"

# --------------------------------------
# arguments for which simulaitons and dates are missing
# --------------------------------------
simulationid=$1 # simulations id, e.g. htm150
from=$2 # from which date the result is missing, e.g. 20241220
to=$3 # to which date the result is missing, e.g. 20241222
# --------------------------------------
# arguments of flexpartset.sh
# --------------------------------------
input=$4 # meteo input path
output=$5 # flexpart output path
id=$6 # species id, NOTE one id for now 
step=$7 # time step in hour
lon=$8 # longitude of release box -180 < LON1 <180
lat=$9 # latitude of release box, -90 < LAT1 < 90
z=${10} # height of release
particles=${11} # Total number of particles to be released

# fix the missing results by each timestep, here is by each day
end=$(date -u -d "$from +1 days" +%Y%m%d)
while [ "$from" != "$to" ]
do
    year=$(date -u -d "$from" +%Y)
    echo simulating $from 

    cd $F_PATH
    # Genereate settings for flexpart simulations
    bash flexpartset_global_14C.sh ${input}${year}/ $output $id $from $end $step $lon $lat $z $particles

    cd $F_PATH/$from 
    # run Flexpart
    nohup $FLEXPARTPATH/FLEXPART pathnames >${output}${from}/log 2>&1 & 

    # Capture the PID of the background process
    PID=$!

    # Wait for the first command to finish
    wait $PID

    cd $F_PATH
    # postprocessing footprint, calculate the delta radiocarbon contribution in permille, and upload result 
    python3 postprocess.py $simulationid $(date -u -d "$from" +%Y-%m-%d)

    from=$(date -u -d "$from +1 days" +%Y%m%d)
    end=$(date -u -d "$end +1 days" +%Y%m%d)
done
exec bash

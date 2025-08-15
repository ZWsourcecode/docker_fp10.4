#!/bin/bash
### run flexpart by date range
### script should be in parent folder of options
### run : bash rundates_eu.sh 
### zhendong.wu@nateko.lu.se

FLEXPARTPATH="/usr/local/flexpart_v10.4/src"
F_PATH="/flexpart"


# --------------------------------------
# arguments for which station and dates
# --------------------------------------
simulationid=$1 # simulations id, e.g. htm150
from=$2 # from which date the result is missing, e.g. 20241220
to=$3 # to which date the result is missing, e.g. 20241222
# --------------------------------------
# arguments of flexpart setting
# --------------------------------------
input=$4 # meteo input path
output=$5 # flexpart output path
id=$6 # species id, NOTE one id for now 
step=$7 # time step in hour
lon=$8 # longitude of release box -180 < LON1 <180
lat=$9 # latitude of release box, -90 < LAT1 < 90
z=${10} # height of release
particles=${11} # Total number of particles to be released

start=${from}
end=${to}

year=$(date -u -d "$from" +%Y)

# Prepare input
bash flexpartset_eu_14Croutine.sh $simulationid ${input}${year}/ $output $id $start $end $step $lon $lat $z $particles
# Change to option directory
cd $F_PATH/${simulationid}/${start} 

# Launch FLEXPART in background and capture PID of the background process
nohup $FLEXPARTPATH/FLEXPART pathnames >${output}${simulationid}/${start}/log 2>&1 & 
PID=$!

echo "FLEXPART started with PID $PID, waiting for it to complete..."
wait $PID
echo "FLEXPART simulation for $start completed."

# Change back to root path
# cd $F_PATH

exec bash

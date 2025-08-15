#!/bin/bash
### automatically run flexpart each day
### script should be in parent folder of options
### run : bash cron_global_14Croutine.sh 9 inputs 
### zhendong.wu@nateko.lu.se

FLEXPARTPATH="/usr/local/flexpart_v10.4/src"
F_PATH="/flexpart"
project="C14"

# ======================================
# arguments 
# ======================================
simulationid=$1 # simulations id, e.g. htm150
input=$2 # meteo input path
# input_nest=${2}${year}/ # meteo nest input path
output=$3 # flexpart output path
id=$4 # species id, NOTE one id for now 
step=$5 # time step in hour
lon=$6 # longitude of release box -180 < LON1 <180
lat=$7 # latitude of release box, -90 < LAT1 < 90
z=$8 # height of release
particles=$9 # Total number of particles to be released

# Remove all digits from input
prefix=${simulationid//[0-9]/}

simulationidproject=${simulationid}${project} # simulations id and project name, e.g. htm150C14eu
output=${output}${prefix}/${simulationidproject}/ # e.g. xxx/htm/htm150C14eu/

date -u

# simulate the day 2 days ago
start=$(date -u -d "-2 days" +%Y%m%d)
end=$(date -u -d "-1 days" +%Y%m%d)
year=$(date -u -d "$start" +%Y)

echo simulating $simulationid $start

cd $F_PATH

# Prepare input and genereate settings for flexpart simulations
bash flexpartset_global_14Croutine.sh $simulationidproject ${input}${year}/ $output $id $start $end $step $lon $lat $z $particles

# Change to option directory
cd $F_PATH/$simulationidproject/$start

# run with mpi
# nohup mpirun -n 2 --allow-run-as-root $FLEXPARTPATH/FLEXPART_MPI pathnames >${output}${start}/log 2>&1 & 

# Launch FLEXPART in background and capture PID of the background process
nohup $FLEXPARTPATH/FLEXPART pathnames >${output}${start}/log 2>&1 & 

# Change back to root path
cd $F_PATH

exec bash


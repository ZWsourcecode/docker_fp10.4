#!/bin/bash
### create pathnames, RELEASES, COMMAND file, and run flexpart
### script should be in parent folder of options
### run : bash run_africa_global.sh 10 inputs 
### zhendong.wu@nateko.lu.se

F_PATH="/flexpart"
# ======================================
# arguments of flexpartset.sh
# ======================================
input=$1 # meteo input path
output=$2 # flexpart output path
id=$3 # species id, NOTE one id for now 
start=$4 # start date, e.g. YYYYMMDD
end=$5 # end date, e.g. YYYYMMDD
step=$6 # time step in hour
lon=$7 # longitude of release box -180 < LON1 <180
lat=$8 # latitude of release box, -90 < LAT1 < 90
z=$9 # height of release
particles=${10} # Total number of particles to be released

FLEXPARTPATH="/usr/local/flexpart_v10.4/src"

# date -u
# start=$(date -u --date="yesterday" +%Y%m%d)
# end=$(date -u +%Y%m%d)

# 5 days later at 14:00 pm using operational data only
# start=$(date -u -d "-2 days" +%Y%m%d)
# end=$(date -u --date="tomorrow" +%Y%m%d)
# end=$(date -u -d "-1 days" +%Y%m%d)
year=$(date -u -d "$start" +%Y)
echo simulating $start 

cd $F_PATH

bash flexpartset_global_africa.sh ${input}${year}/ $output $id $start $end $step $lon $lat $z $particles

cd $F_PATH/$start 
nohup mpirun -n 1 --allow-run-as-root $FLEXPARTPATH/FLEXPART pathnames >${output}${start}/log 2>&1 & 
# nohup $FLEXPARTPATH/FLEXPART pathnames >log 2>&1 & 

cd $F_PATH

# exec bash

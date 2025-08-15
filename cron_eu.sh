#!/bin/bash
### automatically run flexpart each day
### script should be in parent folder of options
### run : bash cron_eu.sh 9 inputs 
### zhendong.wu@nateko.lu.se


# ======================================
# arguments of flexpartset.sh
# ======================================
input=$1 # meteo input path
output=$2 # flexpart output path
id=$3 # species id, NOTE one id for now 
step=$4 # time step in hour
lon=$5 # longitude of release box -180 < LON1 <180
lat=$6 # latitude of release box, -90 < LAT1 < 90
z=$7 # height of release
particles=$8 # Total number of particles to be released

FLEXPARTPATH="/usr/local/flexpart_v10.4/src"
F_PATH="/flexpart"

date -u
# start=$(date -u --date="yesterday" +%Y%m%d)
# end=$(date -u +%Y%m%d)

# 5 days later at 14:00 pm using operational data only
start=$(date -u -d "-3 days" +%Y%m%d)
# end=$(date -u --date="tomorrow" +%Y%m%d)
end=$(date -u -d "-2 days" +%Y%m%d)
year=$(date -u -d "$start" +%Y)
echo simulating $start 

cd $F_PATH

bash flexpartset_eu_14C.sh ${input}${year}/ $output $id $start $end $step $lon $lat $z $particles

cd $F_PATH/$start 
# nohup mpirun -n 1 --allow-run-as-root $FLEXPARTPATH/FLEXPART pathnames >${output}${start}/log 2>&1 & 
nohup $FLEXPARTPATH/FLEXPART pathnames >${output}${start}/log 2>&1 & 


cd $F_PATH

exec bash

#!/bin/bash
### automatically run flexpart each day
### script should be in parent folder of options
### run : bash autorun.sh 9 inputs 
### zhendong.wu@nateko.lu.se

F_PATH="/flexpart"
# ======================================
# arguments of flexpartset.sh
# ======================================
input=$1 # meteo input path
input_nest=$2 # meteo nest input path
output=$3 # flexpart output path
id=$4 # species id, NOTE one id for now 
step=$5 # time step in hour
lon=$6 # longitude of release box -180 < LON1 <180
lat=$7 # latitude of release box, -90 < LAT1 < 90
z=$8 # height of release
particles=$9 # Total number of particles to be released

FLEXPARTPATH="/usr/local/flexpart_v10.4/src"

date -u
# start=$(date -u --date="yesterday" +%Y%m%d)
# end=$(date -u +%Y%m%d)

# two day later at 14:00 pm using operational and forecast data
start=$(date -u -d "-2 days" +%Y%m%d)
# end=$(date -u --date="tomorrow" +%Y%m%d)
end=$(date -u -d "-1 days" +%Y%m%d)
year=$(date -u -d "$start" +%Y)
echo simulating $start 

cd $F_PATH

bash flexpartset_both_14C.sh ${input}${year}/ ${input_nest}${year}/ $output $id $start $end $step $lon $lat $z $particles

cd $F_PATH/$start 
nohup mpirun -n 1 --allow-run-as-root $FLEXPARTPATH/FLEXPART pathnames >${output}${start}/log 2>&1 & 
# nohup $FLEXPARTPATH/FLEXPART pathnames >log 2>&1 & 

cd $F_PATH

exec bash

#!/bin/bash
### automatically run flexpart each day
### script should be in parent folder of options
### run : bash autorun.sh 8 inputs 
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

start="20221121"
end="20221122"
today=$(date -u +%Y%m%d)
# while :; do
while [ "$start" != "$today" ]
do
  echo $start 
  date -u
#   start=$(date -u --date="yesterday" +%Y%m%d)
#   end=$(date -u +%Y%m%d)
  bash flexpartset_global_14C.sh $input $output $id $start $end $step $lon $lat $z $particles
  cd $start
  nohup mpirun -n 1 --allow-run-as-root $FLEXPARTPATH/FLEXPART pathnames >${output}${start}/log 2>&1 & 
  cd ..

  sleep 20m
  start=$(date -u -d "$start +1 days" +%Y%m%d)
  end=$(date -u -d "$end +1 days" +%Y%m%d)
done
exec bash

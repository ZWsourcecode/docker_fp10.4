#!/bin/bash
### automatically run flexpart each day
### script should be in parent folder of options
### run : bash autorun.sh 8 inputs 
### zhendong.wu@nateko.lu.se

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

start="20221120"
end="20221121"
today=$(date -u +%Y%m%d)

# while :; do
while [ "$start" != "$today" ]
do
  echo $start 
  year=$(date -u -d "$start" +%Y)
  date -u
#   start=$(date -u --date="yesterday" +%Y%m%d)
#   end=$(date -u +%Y%m%d)
  bash flexpartset_nest_14C.sh ${input}${year}/ ${input_nest}${year}/ $output $id $start $end $step $lon $lat $z $particles
  cd $start 
  nohup mpirun -n 1 --allow-run-as-root $FLEXPARTPATH/FLEXPART pathnames >${output}${start}/log 2>&1 & 
  cd ..

  sleep 20m
  start=$(date -u -d "$start +1 days" +%Y%m%d)
  end=$(date -u -d "$end +1 days" +%Y%m%d)
done
exec bash

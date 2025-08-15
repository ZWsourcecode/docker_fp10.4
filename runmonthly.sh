#!/bin/bash
### run flexpart monthly
### script should be in parent folder of options
### run : bash runmonthly.sh 
### zhendong.wu@nateko.lu.se

year=2020
declare -a arr_mon=("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12")

# year=2021
# declare -a arr_mon=("01" "02" "03" "04" "05" "06")

# ======================================
# arguments of flexpartset.sh
# ======================================
input=${1}${year}/ # meteo input path
input_nest=${2}${year}/ # meteo nest input path
output=$3 # flexpart output path
id=$4 # species id, NOTE one id for now 
step=$5 # time step in hour
lon=$6 # longitude of release box -180 < LON1 <180
lat=$7 # latitude of release box, -90 < LAT1 < 90
z=$8 # height of release
particles=$9 # Total number of particles to be released


for mon in ${arr_mon[@]}
do
    start=${year}${mon}01
    end=$(date -u -d "$start +1 months" +%Y%m%d)
    
    echo $start 
    date -u

    bash flexpartset_nest_14C.sh $input $input_nest $output $id $start $end $step $lon $lat $z $particles
    cd $start 
    nohup mpirun -n 1 --allow-run-as-root $FLEXPARTPATH/FLEXPART pathnames >${output}${start}/log 2>&1 & 
    cd ..
    
    sleep 15
done
exec bash

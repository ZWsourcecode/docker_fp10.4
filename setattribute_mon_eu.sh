#!/bin/bash
### set global attributes for flexpart footprint
### run : bash setattribute.sh path 
### zhendong.wu@nateko.lu.se

cd $1

for hour in $(ls -d */)
do
    cd $hour

    ncatted -O -h -a backtime,global,o,c,"240 hours" foot
    ncatted -O -h -a unit_conversion,global,o,c,"footprint(s m3 kg-1) / hmix(m) * molar mass of dry air(0.02897 kg mol-1)  = s m2 mol-1 = ppm / (micromoles m-2 s-1)" foot
    ncatted -O -h -a description,global,o,c,"aggregated Flexpart footprints on lon/lat/time grid, aggregated in grid boxes (lat,lon) and Flexpart arriving time (time), aggregated over backtime hours prior to arriving time" foot

    # ncatted -O -h -a backtime,global,o,c,"240 hours" foot_nest
    # ncatted -O -h -a unit_conversion,global,o,c,"footprint(s m3 kg-1) / hmix(m) * molar mass of dry air(0.02897 kg mol-1) = s m2 mol-1 = ppm / (micromoles m-2 s-1)" foot_nest
    # ncatted -O -h -a description,global,o,c,"aggregated Flexpart footprints on lon/lat/time grid, aggregated in grid boxes (lat,lon) and Flexpart arriving time (time), aggregated over backtime hours prior to arriving time" foot_nest

    cd ..
done



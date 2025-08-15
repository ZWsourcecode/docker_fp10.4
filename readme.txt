simulate footprint under mixing height 

The footprint output is stored in NetCDF format, including three main variables: 
1. footprint for each arrival time (spec001_mr_hmix_arr) 
2. accumulated footprints for whole simulation period (spec001_mr_hmix_acc)
3. mixing height (hmix).

Dimension :  
footprint(atime, sfromatime, lon, lat), e.g. the shape is [2881, 240, 180, 360]
hmix(time, lon, lat), e.g. the shape is [3120, 180, 360]
atime: arriving time, e.g. 20180430T230000 to 20171231T230000
sfromatime: seconds from arriving time, e.g. -3600, -7200 … 
time: simulation time, e.g. from 20180430T230000 to 20171222T000000

Footprint: s m3 kg-1
hmix: m


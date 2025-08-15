# Postprocess flexpart footprint
# Convert the unit of footprint from (s m3 kg-1) to (ppm / (micromoles m-2 s-1))
# Organize flexpart footprint by arriving time

import time
from joblib import Parallel, delayed
import xarray as xr
import pandas as pd
# import scipy as sp
import numpy as np
from pathlib import Path
from os import listdir, makedirs, remove, system
from os.path import isfile, join, exists
import sys
import subprocess

# Define function
def grid_area (resolution=0.5):
    """Calculate the area of each grid cell for a user-provided
    grid cell resolution. Area is in square meters, but resolution
    is given in decimal degrees."""
    # Calculations needs to be in radians
    lats = np.deg2rad(np.arange(-90,90+resolution, resolution))
    r_sq = 6371000**2 # % Earth's radius in m^2
    n_lats = int(360./resolution) 
    area = r_sq*np.ones(n_lats)[:, None]*np.deg2rad(resolution)*(
                np.sin(lats[1:]) - np.sin(lats[:-1])) # Spherical cap area, 2 pi r_sq (1-sin(angle))
    return area.T

def parallel_footprint(npoint, IN_PATH, OUT_PATH, darray_fp, darray_hmix, pd_atime, pd_time, ifnest=0):
    """sum footprint during backtime (e.g. 10 days), and save as netcdf per arriving time.
    return 2d dataarray with attributes.
        release point number: npoint , start from 0
        fp unit: s m3 kg-1
        hmix unit: m
        output fp unit: ppm / (micromoles m-2 s-1)
        """
    # slice footprint and hmix for the atime
    fp_10day = darray_fp.isel(atime=npoint)
    # find time index of atime in simulation time, the simulation time and arriveing time may be different 
    hmix_index = pd_time.get_loc(pd_atime[npoint])
    hmix_10day = darray_hmix.isel(time=slice(hmix_index+1,hmix_index+240+1))
    
    try:
        # (s m3 kg-1) / m * (0.02897 kg mol-1) = s m2 mol-1 = ppm / (micromoles m-2 s-1) 
        foot_3d = np.array(fp_10day)/np.array(hmix_10day)*0.02897 
    except:
        foot_3d = np.NaN
    
    # create dataarray
    darry_foot3d = xr.DataArray(
        foot_3d,
        dims=["backtime","lat","lon"],
        coords=dict(
            backtime = pd.to_datetime(np.array(hmix_10day.time)), 
            lat = np.array(hmix_10day.latitude),
            lon = np.array(hmix_10day.longitude)
        ),
        attrs = dict(
            backtime = "240 hours",
            unit_conversion = "footprint(s m3 kg-1) / hmix(m) * molar mass of dry air(0.02897 kg mol-1) = s m2 mol-1 = ppm / (micromoles m-2 s-1)",
            description = "aggregated Flexpart footprints on lon/lat/time grid, aggregated in grid boxes (lat,lon) and Flexpart arriving time (time), aggregated over backtime hours prior to arriving time"
        )
    )
    atime = pd.to_datetime(np.array(fp_10day.atime))
    darry_foot3d = darry_foot3d.assign_coords(time=atime)
    darry_foot3d = darry_foot3d.expand_dims('time')
    
    # sum footprint during backtime
    darry_foot2d = darry_foot3d.sum(dim="backtime")
    
    # add attributes
    darry_foot2d.name = "foot"
    darry_foot2d.attrs["units"] = "ppm per (micromol m-2 s-1)"
    darry_foot2d.attrs["long_name"] = "Flexpart footprints integrated over backward time intervall backtime"

    darry_foot2d.time.attrs['long_name'] = "arriving time"
    darry_foot2d.time.encoding['units'] = f"Hours since {atime.strftime('%Y-%m-%d')} 00:00:00"

    darry_foot2d.lon.attrs['units'] = "degrees_east"
    darry_foot2d.lon.attrs['long_name'] = "longitude in degree east"
    darry_foot2d.lon.attrs['description'] = "grid cell centers"

    darry_foot2d.lat.attrs['units'] = "degrees_north"
    darry_foot2d.lat.attrs['long_name'] = "longitude in degree north"
    darry_foot2d.lat.attrs['description'] = "grid cell centers"
    
    # make directory 
    OUT_FLODER = Station.upper() + "/" + str(atime.year) + "/" + str(atime.month).zfill(2) + "/"
    OUT_FLODER = OUT_FLODER + str(atime.year) + "x" + str(atime.month).zfill(2) + "x" + str(atime.day).zfill(2) + "x" + str(atime.hour).zfill(2)+ "/"
    OUT_FLODER = OUT_PATH + OUT_FLODER

    Path(OUT_FLODER).mkdir(parents=True, exist_ok=True) 
    if ifnest==0:
        file_path_name = OUT_FLODER + "foot"
    else:
        file_path_name = OUT_FLODER + "foot_nest"
    darry_foot2d.to_netcdf(file_path_name)
    
    return OUT_FLODER

def get_footprint(prefix, IN_PATH, OUT_PATH, Station, Project, Year, Month, Day, cpus=1):
    """ get flexpart footprint per arriving time """
    IN_PATH = IN_PATH + Station[0:3] + "/" + Station + Project + "/" + Year + Month + Day + "/"
    
    ds_flexp = xr.open_dataset(IN_PATH+ prefix + '.nc')
    darray_fp = ds_flexp.spec001_mr_hmix_arr
    darray_fp = darray_fp.assign_coords(atime=lambda darray_fp: darray_fp.atime.dt.floor('H'))
    darray_hmix = ds_flexp.hmix
    
    ds_flexp_nest = xr.open_dataset(IN_PATH + prefix + '_nest' + '.nc')
    darray_fp_nest = ds_flexp_nest.spec001_mr_hmix_arr
    darray_fp_nest = darray_fp_nest.assign_coords(atime=lambda darray_fp_nest: darray_fp_nest.atime.dt.floor('H'))
    darray_hmix_nest = ds_flexp_nest.hmix
    
    # atime is arriving time( or particle release time) 
    pd_atime = pd.to_datetime(np.array(darray_fp.atime))
    # time is simulation time
    pd_time = pd.to_datetime(np.array(darray_hmix.time))

# parallel version
    lst_OUT_FLODER = Parallel(n_jobs=cpus,verbose=1)(delayed(parallel_footprint)(n, IN_PATH, OUT_PATH, darray_fp, darray_hmix, pd_atime, pd_time) 
                                     for n in range(len(darray_fp.atime)))
    lst_OUT_FLODER = Parallel(n_jobs=cpus,verbose=1)(delayed(parallel_footprint)(n, IN_PATH, OUT_PATH, darray_fp_nest, darray_hmix_nest, pd_atime, pd_time,ifnest=1) 
                                     for n in range(len(darray_fp.atime)))
    

# ------------------------------------------------
# Extract aggregated Flexpart footprint, and convert unit to  ppm / (micromoles m-2 s-1)
# ------------------------------------------------
IN_PATH = "/data/flexpart/output/"
OUT_PATH = IN_PATH+"flexpartweb/stations/"

Station = sys.argv[1]
date_object = pd.to_datetime(sys.argv[2])

# Extract year, month, and day
Year = date_object.year
Month = date_object.month
Day = date_object.day

Project = "C14"
simulate_date = str(Year)+str(Month)+str(Day).zfill(2)
file_date = str(Year)+str(Month)+str(Day+1).zfill(2)
prefix = "grid_time_" + file_date +"000000"
print("extract footprint of " + Station + " " + str(Year)+str(Month)+str(Day).zfill(2) + " ...")

cpus=1
get_footprint(prefix, IN_PATH, OUT_PATH, Station, Project, str(Year), str(Month), str(Day).zfill(2), cpus=cpus)

ATT_PATH = OUT_PATH + Station.upper() + "/" + str(Year) + "/" + str(Month)

command = ["./setattribute_mon.sh", ATT_PATH]

system("chmod +x setattribute_mon.sh")
# Run the command
try:
    subprocess.run(command, check=True)
except subprocess.CalledProcessError as e:
    print(f"Error: {e}")

# ------------------------------------------------
# Calculate delta radiocarbon contribution, permil
# ------------------------------------------------
print("Calculate delta radiocarbon of " + Station + " " + str(Year)+str(Month)+str(Day).zfill(2) + " ...")

# setting parameter for calculating delta radiocarbon contribution
url_noaa_co2_mm_gl = 'https://gml.noaa.gov/webdata/ccgg/trends/co2/co2_mm_gl.txt'
co2_global=pd.read_csv(url_noaa_co2_mm_gl,header=None,
                     names=['year','month', 'decimal', 'average','average_unc','trend','trend_unc'], skiprows=40,delimiter=r"\s+")

Xco2 = co2_global['average'].iloc[-1] # NOAA global monthly Average CO2, ppm
Mc = 12 # g/mol
Aabs = 0.226 # Bq/gC

# Radiocarbon (14CO2) emissions from nuclear facilities in 2020
# The data is derived from the European Commission RAdioactive Discharges Database(RADD, Zazzeri et al. (2018)).
PATH_NUCLEAR = IN_PATH+"flux/"
darray_Q = xr.open_dataarray(PATH_NUCLEAR+"Radiocarbon_nuclear_emissions_2021_eu.nc")
darray_Q = darray_Q.astype(np.float64)

lst_QF = [] 
delta_14C = pd.DataFrame()
delta_14C["UTC"] = []
delta_14C["14C"] = []
for Hour in range(24):
    floder_Flexpart = str(Year)+"x"+str(Month)+"x"+str(Day).zfill(2)+"x"+str(Hour).zfill(2)
    darray_flexpart = xr.open_dataarray(ATT_PATH+"/"+floder_Flexpart+"/foot_nest", engine='netcdf4')
    darray_flexpart = darray_flexpart.astype(np.float64)
    # QF: ppm Bq micromol-1
    darray_QF = np.array(darray_Q)*np.array(darray_flexpart)[0,:,:]
    # QF: ppm Bq mol-1
    QF = darray_QF.sum()*1e6
    lst_QF.append(QF)
    delta_14C.loc[Hour,"UTC"] = pd.to_datetime(np.array(darray_flexpart.time))
    delta_14C.loc[Hour,"14C"] = round(QF/(Xco2 * Mc * Aabs) * 1000,3)
delta_14C["ifkeep"] = delta_14C["14C"] < 0.5
delta_14C.to_csv(ATT_PATH+"/" +"delta_14C_" + Station + "_" + simulate_date + ".csv", header=True,index=False, na_rep= "NaN")
delta_14C
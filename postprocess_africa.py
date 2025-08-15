# Postprocess flexpart footprint
# Convert the unit of footprint from (s m3 kg-1) to (ppm / (micromoles m-2 s-1))
# Organize flexpart footprint by arriving time
# Calculate delta radiocarbon contribution, permil

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
# import pysftp
import paramiko
import requests

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

def loop_footprint (npoint, IN_PATH, OUT_PATH, darray_fp, darray_hmix, pd_atime, pd_time, ifnest=0):
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
    # lst_OUT_FLODER = Parallel(n_jobs=cpus,verbose=1)(delayed(parallel_footprint)(n, IN_PATH, OUT_PATH, darray_fp, darray_hmix, pd_atime, pd_time) 
    #                                  for n in range(len(darray_fp.atime)))
    # lst_OUT_FLODER = Parallel(n_jobs=cpus,verbose=1)(delayed(parallel_footprint)(n, IN_PATH, OUT_PATH, darray_fp_nest, darray_hmix_nest, pd_atime, pd_time,ifnest=1) 
    #                                  for n in range(len(darray_fp.atime)))
    
    # normal loop version    
    for npoint in range(len(darray_fp.atime)):
        loop_footprint(npoint, IN_PATH, OUT_PATH, darray_fp, darray_hmix, pd_atime, pd_time)
        loop_footprint(npoint, IN_PATH, OUT_PATH, darray_fp_nest, darray_hmix_nest, pd_atime, pd_time, ifnest=1)
        if npoint % 10 ==0:
            print(npoint, end="...")
    print("done")

def upload_to_sftp(hostname, port, username, password, local_file_path, remote_directory, remote_file_name):
    # Create an SSH client
    ssh = paramiko.SSHClient()

    # Automatically add the server's host key (this is insecure; see comments below)
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())

    try:
        # Connect to the SFTP server
        ssh.connect(hostname, port, username, password)

        # Create an SFTP client
        sftp = ssh.open_sftp()

        try:
            # Create the remote directory if it doesn't exist
            current_path = '/'
            # Create the remote directory if it doesn't exist
            dirs = remote_directory.split('/')
            for dir_name in dirs:
                if dir_name:
                    current_path = current_path + dir_name + '/'
                    try:
                        sftp.stat(current_path)
                    except IOError as e:
                        sftp.mkdir(current_path)
                        print(f"Created directory: {current_path}")
            # try:
            #     sftp.stat(remote_directory)
            # except FileNotFoundError:
            #     sftp.mkdir(remote_directory)

            # Upload the local file to the remote directory
            sftp.put(local_file_path, f"{remote_directory}/{remote_file_name}")
                # Change to the remote directory
            # sftp.cwd(remote_directory)

            # # Upload the file
            # sftp.put(local_file_path, remote_file_name)

            print(f"File '{local_file_path}' uploaded to '{remote_directory}' on {hostname}")

        finally:
            # Close the SFTP connection
            sftp.close()

    finally:
        # Close the SSH connection
        ssh.close()    

# ------------------------------------------------
# Extract aggregated Flexpart footprint, and convert unit to  ppm / (micromoles m-2 s-1)
# ------------------------------------------------
IN_PATH = "/flexpart/postprocess/"
OUT_PATH = IN_PATH+"flexpartweb/Africa/"

Station = sys.argv[1]
date_object = pd.to_datetime(sys.argv[2])

# Extract year, month, and day
Year = str(date_object.year)
Month = str(date_object.month).zfill(2)
Day = str(date_object.day).zfill(2)

Project = ""

ifrunbymon = 1 # 1: simulation run by month, 0: simulation run by day
cpus=1

pattern = Year+Month
datelist = [fn for fn in listdir(IN_PATH+Station[0:3] + "/" + Station + Project) if pattern in fn]
datelist.sort()
for date in datelist:
    yy = date[0:4]
    mm = date[4:6]
    dd = date[6:8]
    simulate_date = date
    if ifrunbymon == 0:  
        file_date = pd.to_datetime(date) + pd.DateOffset(days=1)
    else:
        file_date = pd.to_datetime(date) + pd.DateOffset(months=1)
    file_date = file_date.strftime('%Y%m%d')
    prefix = "grid_time_" + file_date +"000000"
    print("extract footprint for " + Station + Project + " " + simulate_date + " ...")
    get_footprint(prefix, IN_PATH, OUT_PATH, Station, Project, yy, mm, dd, cpus=cpus)

ATT_PATH = OUT_PATH + Station + "/" + Year + "/" + Month
command = ["/flexpart/setattribute_mon.sh", ATT_PATH]

# system("chmod +x setattribute_mon.sh")
# Run the command
try:
    subprocess.run(command, check=True)
except subprocess.CalledProcessError as e:
    print(f"Error: {e}")

print("Done")

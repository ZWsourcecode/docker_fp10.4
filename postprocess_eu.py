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
    
    # ds_flexp_nest = xr.open_dataset(IN_PATH + prefix + '_nest' + '.nc')
    # darray_fp_nest = ds_flexp_nest.spec001_mr_hmix_arr
    # darray_fp_nest = darray_fp_nest.assign_coords(atime=lambda darray_fp_nest: darray_fp_nest.atime.dt.floor('H'))
    # darray_hmix_nest = ds_flexp_nest.hmix
    
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
        # loop_footprint(npoint, IN_PATH, OUT_PATH, darray_fp_nest, darray_hmix_nest, pd_atime, pd_time, ifnest=1)
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
OUT_PATH = IN_PATH+"flexpartweb/stations_eu/"

Station = sys.argv[1]
date_object = pd.to_datetime(sys.argv[2])

# Extract year, month, and day
Year = date_object.year
Month = date_object.month
Day = date_object.day

Project = "C14eu"
simulate_date = str(Year)+str(Month).zfill(2)+str(Day).zfill(2)
# file_date = str(Year)+str(Month).zfill(2)+str(Day+1).zfill(2)

date_object = pd.to_datetime(simulate_date, format="%Y%m%d")
new_date_object = date_object + pd.Timedelta(days=1)
file_date = new_date_object.strftime("%Y%m%d")
file_day = pd.to_datetime(file_date, format="%Y%m%d").day

prefix = "grid_time_" + file_date +"000000"
print("extract footprint of " + Station + " " + str(Year)+str(Month).zfill(2)+str(Day).zfill(2) + " ...")

cpus=1
get_footprint(prefix, IN_PATH, OUT_PATH, Station, Project, str(Year), str(Month).zfill(2), str(Day).zfill(2), cpus=cpus)

ATT_PATH = OUT_PATH + Station.upper() + "/" + str(Year) + "/" + str(Month).zfill(2)

command = ["/flexpart/setattribute_mon_eu.sh", ATT_PATH]

# system("chmod +x setattribute_mon.sh")
# Run the command
try:
    subprocess.run(command, check=True)
except subprocess.CalledProcessError as e:
    print(f"Error: {e}")

print("Done")

# ------------------------------------------------
# Calculate delta radiocarbon contribution, permil
# ------------------------------------------------
print("Calculate delta radiocarbon of " + Station + " " + str(Year)+str(Month).zfill(2)+str(Day).zfill(2) + " ...")

# setting parameter for calculating delta radiocarbon contribution
url_noaa_co2_mm_gl = 'https://gml.noaa.gov/webdata/ccgg/trends/co2/co2_mm_gl.txt'

# On the first day of the month, download/update the co2 file which will be used in case the NOAA server is down.
if file_day == 1:
    # Send a HTTP GET request to the URL
    response = requests.get(url_noaa_co2_mm_gl)

    # Check if the request was successful
    if response.status_code == 200:
        # Save the content to a local file
        with open(IN_PATH + "co2_mm_gl.txt", "wb") as file:
            file.write(response.content)
        print("File downloaded and saved as co2_mm_gl.txt")
    else:
        print(f"Failed to download the file. Status code: {response.status_code}")

# try to read the noaa co2 trend online, if NOAA server is down or can not read the data, use the local data instead
try:
    co2_global=pd.read_csv(url_noaa_co2_mm_gl,header=None,
                     names=['year','month', 'decimal', 'average','average_unc','trend','trend_unc'], skiprows=40,delimiter=r"\s+")
except:
    print("Can not read the co2 file remotely, use local file instead")
    co2_global=pd.read_csv(IN_PATH + "co2_mm_gl.txt",header=None,
                     names=['year','month', 'decimal', 'average','average_unc','trend','trend_unc'], skiprows=40,delimiter=r"\s+")

Xco2 = co2_global['average'].iloc[-1] # NOAA global monthly Average CO2, ppm
Mc = 12 # g/mol
Aabs = 0.226 # Bq/gC

# Radiocarbon (14CO2) emissions from nuclear facilities in 2020
# The data is derived from the European Commission RAdioactive Discharges Database(RADD, Zazzeri et al. (2018)).

# PATH_NUCLEAR = IN_PATH+"flux/"
# darray_Q = xr.open_dataarray(PATH_NUCLEAR+"Radiocarbon_nuclear_emissions_2022_eu.nc")
# darray_Q = darray_Q.astype(np.float64)

# lst_QF = [] 
# delta_14C = pd.DataFrame()
# delta_14C["UTC"] = []
# delta_14C["14C"] = []
# for Hour in range(24):
#     floder_Flexpart = str(Year)+"x"+str(Month).zfill(2)+"x"+str(Day).zfill(2)+"x"+str(Hour).zfill(2)
#     darray_flexpart = xr.open_dataarray(ATT_PATH+"/"+floder_Flexpart+"/foot", engine='netcdf4')
#     darray_flexpart = darray_flexpart.astype(np.float64)
#     # QF: ppm Bq micromol-1
#     darray_QF = np.array(darray_Q)*np.array(darray_flexpart)[0,:,:]
#     # QF: ppm Bq mol-1
#     QF = darray_QF.sum()*1e6
#     lst_QF.append(QF)
#     delta_14C.loc[Hour,"UTC"] = np.datetime_as_string(darray_flexpart.time, unit='s').item()
#     delta_14C.loc[Hour,"14C"] = round(QF/(Xco2 * Mc * Aabs) * 1000,3)
# delta_14C["ifkeep"] = delta_14C["14C"] < 0.5
# Filename = "delta_14C_" + Station + "_" + simulate_date + "eu.csv"
# delta_14C.to_csv(ATT_PATH+"/" + Filename, header=True,index=False, na_rep= "NaN")

# unclear emission based on the year
strYear = '2022'
PATH_NUCLEAR = IN_PATH+"flux/"
df_nuclear = pd.read_csv(PATH_NUCLEAR+"emission_dataset_2023_10_27.csv")
df_nuclear.replace(',', '.', regex=True, inplace=True)
df_nuclear[strYear] = pd.to_numeric(df_nuclear[strYear], errors='coerce')
df_nuclear['lat'] = pd.to_numeric(df_nuclear['lat'], errors='coerce')
df_nuclear['lon'] = pd.to_numeric(df_nuclear['lon'], errors='coerce')

lon_global = np.arange(-179.95,180,0.1)
lat_global = np.arange(-89.95,90,0.1)
da_global_14co2 = xr.DataArray(0.0, coords=[lat_global,lon_global], dims= ["lat","lon"])
da_eu_14co2 = da_global_14co2.sel(lon=slice(-15,35),lat=slice(33,73))

da_global_area = xr.DataArray(grid_area(0.1), coords=[lat_global,lon_global], dims= ["lat","lon"])
da_eu_area = da_global_area.sel(lon=slice(-15,35),lat=slice(33,73))

Xco2 = co2_global['average'].iloc[-1] # NOAA global monthly Average CO2, ppm
Mc = 12 # g/mol
Aabs = 0.226 # Bq/gC

df_foot = df_nuclear[['lat','lon','country','facility',strYear]].copy()
df_delta14C_NPP = df_nuclear[['lat','lon','country','facility',strYear]].copy()
delta_14C = pd.DataFrame()
delta_14C["datetime"] = []
delta_14C["datetime"] = delta_14C["datetime"].astype(str)
delta_14C["14C"] = []
Hour_columns = [str(i).zfill(2) for i in range(24)] 
for strHour in Hour_columns:
    df_foot[strHour] = 0.0
    df_delta14C_NPP[strHour] = 0.0
    floder_Flexpart = str(Year)+"x"+str(Month).zfill(2)+"x"+str(Day).zfill(2)+"x"+strHour
    darray_flexpart = xr.open_dataarray(ATT_PATH+"/"+floder_Flexpart+"/foot")
    darray_flexpart = darray_flexpart.astype(np.float64)
    delta_14C.loc[strHour,"datetime"] = np.datetime_as_string(darray_flexpart.time, unit='s').item()
    delta_14C.loc[strHour,"14C"] = 0.0
    for i in range(df_foot.shape[0]):
        lon = df_foot.lon[i]
        lat = df_foot.lat[i]
        # footprint, ppm per (micromol m-2 s-1)
        F = darray_flexpart.sel(lon=lon, lat=lat, method='nearest')
        df_foot.loc[i,strHour] = F
        # grid area, m2
        area = da_eu_area.sel(lon=lon, lat=lat, method='nearest')
        # GBq yr-1 to Bq m-2 s-1, GBq x 1e9 /(365 x 24 x 60 x 60) / area(m2) = Bq m-2 s-1
        # Q, Bq m-2 s-1
        # Q = darray_Q.sel(lon=lon, lat=lat, method='nearest')
        Q = df_delta14C_NPP.loc[i,strYear]*1e9/(365*24*60*60)/area.values
        # QF: ppm Bq mol-1
        QF = Q*F*1e6
        QF = QF.values.item()
        # delta14C: permil
        df_delta14C_NPP.loc[i,strHour] = round(QF/(Xco2 * Mc * Aabs) * 1000,5)
        delta_14C.loc[strHour,"14C"] = delta_14C.loc[strHour,"14C"] + round(QF/(Xco2 * Mc * Aabs) * 1000,5)
delta_14C["ifkeep"] = delta_14C["14C"] < 0.5
simulate_date = str(Year)+str(Month).zfill(2)+str(Day).zfill(2)
Filename = "delta_14C_" + Station + "_" + simulate_date + "eu.csv"
delta_14C.to_csv(ATT_PATH+"/" + Filename, header=True,index=False, na_rep= "NaN")
# Filename = "delta_14C_" + Station + "_" + simulate_date + "eu_NPP.csv"
# df_delta14C_NPP.to_csv(ATT_PATH+"/" + Filename, header=True,index=False, na_rep= "NaN")

print("Calculation of delta 14C is done")


# ------------------------------------------------
# upload data
# ------------------------------------------------
print("Start uploading ...")
# SFTP Connection Information
host = 'icos-ssh.lsce.ipsl.fr'
username = 'sampling'
password = 'zBTp6582y'
port = 22  # default SFTP port

# Local file path
local_file_path = ATT_PATH+"/" + Filename

# Remote directory and file name
remote_directory = '/data/CPrequests/icoscp_eu/'+Station +"/" + str(Year)+"/"+str(Month).zfill(2)
remote_file_name = Filename


# # Connect to the SFTP server
# with pysftp.Connection(host, username=username, password=password, port=port) as sftp:
#     print(f"Connected to {host}")

#     # Check if the remote directory exists; if not, create it
#     if not sftp.exists(remote_directory):
#         sftp.makedirs(remote_directory)
#         print(f"Created remote directory: {remote_directory}")

#     # Change to the remote directory
#     sftp.cwd(remote_directory)

#     # Upload the file
#     sftp.put(local_file_path, remote_file_name)
#     print(f"File '{local_file_path}' uploaded to '{remote_directory}/{remote_file_name}'")

# print("SFTP connection closed.")

upload_to_sftp(host, port, username, password, local_file_path, remote_directory,remote_file_name)

print("Done")

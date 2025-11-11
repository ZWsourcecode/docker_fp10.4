import glob
import os
import pandas as pd
import datetime as dt
from pathlib import Path
import numpy as np

def merge_month(in_path, station, date, domain="global"):
    """Merge daily delta_14C CSV files for a station/month into a single monthly CSV.
    Args:
        in_path: Path or string to data root (e.g. app_dir / 'data/v11')
        station: station code (case-insensitive)
        date: a datetime.date or datetime-like with .year and .month
        domain: 'global' or 'eu' to select file name pattern
    Returns:
        pd.DataFrame: concatenated and sorted DataFrame (empty DataFrame if no files).
    """

    in_path = Path(in_path)
    station_str = str(station)
    station_lower = station_str.lower()
    Year = str(date.year)
    Month = str(date.month).zfill(2)

    # Build pattern and outfile depending on domain
    if domain == "global":
        datapath = in_path / "stations" / station_str / Year / Month
        outpath = in_path / "stations" / station_str / Year
        pattern = f"{datapath}/**/delta_14C_{station_lower}_[0-9]*.csv"
        outfile = outpath / f"delta_14C_{station_lower}_month_{Year}{Month}.csv"
    elif domain == "eu":
        datapath = in_path / "stations_eu" / station_str / Year / Month
        outpath = in_path / "stations_eu" / station_str / Year
        pattern = f"{datapath}/**/delta_14C_{station_lower}_[0-9]*eu.csv"
        outfile = outpath / f"delta_14C_{station_lower}_month_{Year}{Month}eu.csv"
    else:
        raise ValueError("domain must be 'global' or 'eu'")
    outpath.mkdir(parents=True, exist_ok=True)

    # Find all matching CSV files recursively
    csv_files = glob.glob(pattern, recursive=True)

    if not csv_files:
        print(f"No files found for station={station_str}, path={datapath}")
        return pd.DataFrame()

    all_dataframes = []
    for csv_file in sorted(csv_files):
        try:
            df = pd.read_csv(csv_file)
            # Add source file info if needed
            df['source_file'] = os.path.basename(csv_file)
            all_dataframes.append(df)
        except Exception as e:
            print(f"Error reading {csv_file}: {e}")
            continue

    if not all_dataframes:
        print("No readable CSV files found.")
        return pd.DataFrame()

    # Concatenate all DataFrames
    combined_df = pd.concat(all_dataframes, ignore_index=True)

    # Normalize and sort by datetime if present
    col = 'datetime'
    if col in combined_df.columns:
        try:
            combined_df[col] = pd.to_datetime(combined_df[col], errors='coerce')
            combined_df = combined_df.sort_values(col).reset_index(drop=True)
        except Exception as e:
            print(f"Error processing {col}: {e}")

    # Save the combined DataFrame to a new CSV file
    try:
        combined_df.to_csv(outfile, index=False)
    except Exception as e:
        print(f"Could not write output file {outfile}: {e}")
    print(f"Merged {len(csv_files)} files into {os.path.basename(outfile)}")
    # return combined_df

def merge_months(in_path, station, date, domain="global"):
    """For each month in the year, merge daily delta_14C CSV files for a station/month into a single monthly CSV.
    Args:
        in_path: Path or string to data root (e.g. app_dir / 'data/v11')
        station: station code (case-insensitive)
        date: a datetime.date or datetime-like with .year used to select the year
        domain: 'global' or 'eu' to select which station subfolder to check ('stations' vs 'stations_eu')
    """

    in_path = Path(in_path)
    station_str = str(station)
    Year = str(date.year)

    # Choose base folder depending on domain to match merge_month behavior
    base_folder = "stations" if domain == "global" else "stations_eu"
    for month in range(1, 13):
        month_z = str(month).zfill(2)
        monthpath = in_path / base_folder / station_str / Year / month_z
        # print(f"{monthpath}, monthpath.exists()={monthpath.exists()}")
        if not monthpath.exists():
            continue  # Skip months with no data directory
        month_date = date.replace(month=month, day=1)
        merge_month(in_path, station, month_date, domain=domain)

def merge_year(in_path, station, date, domain="global", overwrite=False):
    """
    Merge the daily CSVs (or create them) in a year into a single yearly CSV.
    Args:
        in_path: Path or string to data root (e.g. app_dir / 'data/v11')
        station: station code (case-insensitive)
        date: a datetime.date or datetime-like with .year and .month
        domain: 'global' or 'eu' to select file name pattern
    Returns:
        pd.DataFrame: concatenated and sorted DataFrame (empty DataFrame if no files).
    """

    in_path = Path(in_path)
    station_str = str(station)
    station_lower = station_str.lower()
    Year = str(date.year)

    # Determine yearly outfile name based on domain
    if domain == "global":
        datapath = in_path / "stations" / station_str / Year
        outpath = in_path / "stations" / station_str / Year
        pattern = f"{datapath}/**/delta_14C_{station_lower}_[0-9]*.csv"
        outfile = outpath / f"delta_14C_{station_lower}_year_{Year}.csv"
    elif domain == "eu":
        datapath = in_path / "stations_eu" / station_str / Year 
        outpath = in_path / "stations_eu" / station_str / Year
        pattern = f"{datapath}/**/delta_14C_{station_lower}_[0-9]*eu.csv"
        outfile = outpath / f"delta_14C_{station_lower}_year_{Year}eu.csv"
    else:
        raise ValueError("domain must be 'global' or 'eu'")
    outpath.mkdir(parents=True, exist_ok=True)

    if outfile.exists() and not overwrite:
        # Load and return existing file to avoid recomputing
        try:
            return pd.read_csv(outfile)
        except Exception as e:
            print(f"Could not read existing yearly file {outfile}: {e}")
            # continue and attempt to rebuild
    
    # Find all matching CSV files recursively
    csv_files = glob.glob(pattern, recursive=True)

    if not csv_files:
        print(f"No files found for station={station_str}, path={datapath}")
        return pd.DataFrame()

    all_dataframes = []
    for csv_file in sorted(csv_files):
        try:
            df = pd.read_csv(csv_file)
            # Add source file info if needed
            df['source_file'] = os.path.basename(csv_file)
            all_dataframes.append(df)
        except Exception as e:
            print(f"Error reading {csv_file}: {e}")
            continue

    if not all_dataframes:
        print("No readable CSV files found.")
        return pd.DataFrame()

    # Concatenate all DataFrames
    combined_df = pd.concat(all_dataframes, ignore_index=True)

    # Sort by datetime if present
    col = 'datetime'
    if col in combined_df.columns:
        try:
            combined_df[col] = pd.to_datetime(combined_df[col], errors='coerce')
            combined_df = combined_df.sort_values(col).reset_index(drop=True)
        except Exception as e:
            print(f"Error processing {col}: {e}")

    # Save the combined DataFrame to a new CSV file
    try:
        combined_df.to_csv(outfile, index=False)
    except Exception as e:
        print(f"Could not write output file {outfile}: {e}")
    print(f"Merged {len(csv_files)} files into {os.path.basename(outfile)}")
    # return combined_df


#!/usr/bin/env python3
"""
Script to generate monthly 14C data by merging data for all stations.
Processes both global and EU domain stations.

Usage:
    python gen_month14C.py [--date YYYY-MM-DD] [--domain global|eu|both] [--station STATION_CODE]
    
Examples:
    python gen_month14C.py --date 2025-08-01 --domain global
    python gen_month14C.py --station hfd100 --domain eu
"""
from shared import merge_month, merge_months, merge_year
import datetime as dt
import pandas as pd
import os
import argparse
import sys
from pathlib import Path

# Configuration
IN_PATH = "/flexpart/postprocess/"
MONTH_PATH = IN_PATH + "flexpartweb/"
DEFAULT_DATE = "2025-08-01"

# Station configuration files
STATION_CONFIGS = {
    "global": "stations.conf",
    "eu": "stations_eu.conf"
}


def load_station_config(config_file):
    """
    Load station configuration from a text file with space-separated values.
    
    Args:
        config_file (str): Path to the station configuration file
        
    Returns:
        pd.DataFrame: DataFrame containing station information
    """
    if not os.path.exists(config_file):
        raise FileNotFoundError(f"Station config file not found: {config_file}")
    
    with open(config_file, "r") as f:
        lines = [line.strip() for line in f if line.strip() and not line.strip().startswith('#')]
    
    # Parse the space-separated values into a DataFrame
    data = []
    for line in lines:
        parts = line.split()
        if len(parts) >= 5:  # Ensure we have all required columns
            data.append(parts)
        else:
            print(f"Warning: Skipping malformed line: {line}", flush=True)
    
    if not data:
        raise ValueError(f"No valid station data found in {config_file}")
    
    return pd.DataFrame(data, columns=["station", "lon", "lat", "z", "cpuid"])


def process_stations(config_file, domain, target_date, specific_station=None):
    """
    Process all stations for a given domain.
    
    Args:
        config_file (str): Path to the station configuration file
        domain (str): Domain name ('global' or 'eu')
        target_date (datetime.date): Date for processing
        specific_station (str, optional): Process only this specific station
    """
    print(f"\n=== Processing {domain.upper()} domain stations ===", flush=True)
    
    try:
        station_conf = load_station_config(config_file)
        stations = station_conf["station"].tolist()
        
        # Filter to specific station if requested
        if specific_station:
            stations = [s for s in stations if s.lower() == specific_station.lower()]
            if not stations:
                print(f"Station '{specific_station}' not found in {config_file}", flush=True)
                return
        
        print(f"Found {len(stations)} station(s) to process in {config_file}", flush=True)
        
        success_count = 0
        error_count = 0
        
        for i, station in enumerate(stations, 1):
            print(f"Processing station {i}/{len(stations)}: {station}", flush=True)
            try:
                merge_months(MONTH_PATH, station.upper(), target_date, domain=domain)
                print(f"  ✓ Successfully processed {station}", flush=True)
                success_count += 1
            except Exception as e:
                print(f"  ✗ Error processing {station}: {e}", flush=True)
                error_count += 1
        
        print(f"\n{domain.upper()} Summary: {success_count} successful, {error_count} errors", flush=True)
                
    except Exception as e:
        print(f"Error processing {domain} stations: {e}", flush=True)


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate monthly 14C data by merging data for stations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --date 2025-08-01 --domain global
  %(prog)s --station hfd100 --domain eu
  %(prog)s --date 2024-12-01 --domain both
        """
    )
    
    parser.add_argument(
        '--date', '-d',
        type=str,
        default=DEFAULT_DATE,
        help=f'Date in YYYY-MM-DD format (default: {DEFAULT_DATE})'
    )
    
    parser.add_argument(
        '--domain',
        choices=['global', 'eu', 'both'],
        default='both',
        help='Domain to process: global, eu, or both (default: both)'
    )
    
    parser.add_argument(
        '--station', '-s',
        type=str,
        help='Process only a specific station code (case-insensitive)'
    )
    
    parser.add_argument(
        '--list-stations',
        action='store_true',
        help='List all available stations and exit'
    )
    
    return parser.parse_args()


def list_stations():
    """List all available stations from both config files."""
    print("Available stations:", flush=True)
    
    for domain, config_file in STATION_CONFIGS.items():
        try:
            station_conf = load_station_config(config_file)
            stations = station_conf["station"].tolist()
            print(f"\n{domain.upper()} domain ({len(stations)} stations):", flush=True)
            for i, station in enumerate(stations, 1):
                print(f"  {i:2d}. {station}", flush=True)
        except Exception as e:
            print(f"\n{domain.upper()} domain: Error loading {config_file} - {e}", flush=True)


def main():
    """Main function to process station configurations."""
    args = parse_arguments()
    
    # List stations and exit if requested
    if args.list_stations:
        list_stations()
        return
    
    # Parse and validate the date
    try:
        target_date = dt.datetime.strptime(args.date, "%Y-%m-%d").date()
    except ValueError as e:
        print(f"Error: Invalid date format '{args.date}'. Use YYYY-MM-DD format.", flush=True)
        return
    
    print(f"Processing data for date: {target_date}", flush=True)
    if args.station:
        print(f"Processing only station: {args.station}", flush=True)
    
    # Determine which domains to process
    domains_to_process = []
    if args.domain == 'both':
        domains_to_process = ['global', 'eu']
    else:
        domains_to_process = [args.domain]
    
    # Process each domain
    total_success = 0
    total_errors = 0
    
    for domain in domains_to_process:
        config_file = STATION_CONFIGS[domain]
        try:
            process_stations(config_file, domain, target_date, args.station)
        except Exception as e:
            print(f"Critical error processing {domain} domain: {e}", flush=True)
    
    print(f"\n=== Processing completed for {target_date} ===", flush=True)
    sys.stdout.flush()  # Ensure all output is written


if __name__ == "__main__":
    main()
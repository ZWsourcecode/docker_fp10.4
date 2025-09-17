#!/bin/bash

# --- Configuration ---
INPUT="/meteo/eur/"
OUTPUT="/flexpart/output/"
SCRIPT="/flexpart/cron_eu_14Croutine.sh"
POSTSCRIPT="/flexpart/postprocess_eu.py"

ID=24              # Species ID (currently only supports one)
STEP=1             # Time step in hours
PARTICLES=200      # Number of particles to be released

START_HOUR=0      # First job starts at 0:05
START_MIN=5

# --- Initial Time Setup ---
MINUTE=$START_MIN
HOUR=$START_HOUR

# --- Process Each Station ---
# stations.conf format: station_name lon lat z cpuid
# Example line: hfd100 0.231 50.977 100 2

while read -r station lon lat z cpuid; do

    # Skip comments and empty lines
    [[ "$station" =~ ^#.*$ || -z "$station" ]] && continue  

    # Simulation crontab line
    echo "$MINUTE $HOUR * * * taskset -c $cpuid /bin/bash $SCRIPT $station $INPUT $OUTPUT $ID $STEP $lon $lat $z $PARTICLES >> /flexpart/logcron$station"

    # Postprocessing 1 hour later
    POST_HOUR=$((HOUR + 1))
    (( POST_HOUR >= 24 )) && POST_HOUR=0

    echo "$MINUTE $POST_HOUR * * * taskset -c $cpuid /usr/bin/python3 $POSTSCRIPT $station \$(date -u -d \"-3 days\" +\\%Y-\\%m-\\%d) >> /flexpart/logpython$station"
    echo

    # Increment time by 10 minutes
    MINUTE=$((MINUTE + 10))
    if (( MINUTE >= 60 )); then
        MINUTE=$((MINUTE - 60))
        HOUR=$((HOUR + 1))
        (( HOUR >= 24 )) && HOUR=0
    fi

done < stations_eu.conf
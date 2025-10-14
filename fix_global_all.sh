#!/bin/bash
### run : bash fix_global_all.sh START_DATE END_DATE 

# Define the fixed parameters
START_DATE=$1
END_DATE=$2
INPUT_DIR="/meteo/"
OUTPUT_DIR="/flexpart/postprocess/"
SCRIPT="/flexpart/fix_global_14Croutine.sh"
F_PATH="/flexpart"

ID=24              # Species ID (currently only supports one)
STEP=1             # Time step in hours
PARTICLES=500      # Number of particles to be released

# Define arguments for each station
ARGS_LIST=(
  "hfd100 0.231 50.977 100 1"
  "htm150 13.4189 56.0976 150 2"
  "gat341 11.4429 53.0657 341 3"
  "kre250 15.08000 49.57200 250 4"
  "lin098 14.12260 52.16630 98 5"
  "ope120 5.50360 48.56190 120 6"
  "ste252 8.45880 53.04310 252 7"
  "trn180 2.11250 47.96470 180 8"
  "hpb131 11.02460 47.80110 131 9"
  "kit200 8.42490 49.09150 200 10"
  "cbw200 4.926161 51.970258 200 11"
  "bik300 23.026750 53.231528 300 12"
  "jfj960 7.9851 46.5475 960 13"
  "oxk163 11.80830 50.03000 163 14"
  "ssl120 7.9113 47.9106 120 15"
  "pal012 24.11570 67.97330 12 16"
  "sac100 2.14200 48.72270 100 17"
  "lmp008 12.63220 35.51810 8 18"
  "pot100 15.7237 40.601 100 19"
  "zep015 11.88670 78.90720 15 20"
  "nor100 17.47940 60.08640 100 21"
  "svb150 19.77500 64.25600 150 22"
)


# Loop through and process each station
for ARGS in "${ARGS_LIST[@]}"; do
    # Parse arguments
    read -r station lon lat z cpuid <<< "$ARGS"

    # Skip comments and empty lines
    # [[ "$station" =~ ^#.*$ || -z "$station" ]] && continue  

    # Debugging output
    echo "Processing station: $station, lon: $lon, lat: $lat, z: $z, cpuid: $cpuid"

    # Check if all required fields are present
    if [[ -z "$station" || -z "$lon" || -z "$lat" || -z "$z" || -z "$cpuid" ]]; then
        echo "Error: Missing fields in stations.conf for station $station. Skipping..."
        continue
    fi

    echo "Launching fix missing from $START_DATE to $END_DATE for $station on CPU core $cpu..."
    
    # Run the job in the background using taskset + nohup
    # taskset -c $cpuid nohup bash $SCRIPT $station $START_DATE $END_DATE $INPUT_DIR $OUTPUT_DIR $ID $STEP $lon $lat $z $PARTICLES > logfix_${station} 2>&1 &

    # Run the job in the background using taskset + nohup
    nohup bash -c "taskset -c $cpuid bash $SCRIPT $station $START_DATE $END_DATE $INPUT_DIR $OUTPUT_DIR $ID $STEP $lon $lat $z $PARTICLES" > logfix_${station} 2>&1 &

    # Wait 10 minutes before starting the next one
    echo "Sleeping for 10 minutes before launching the next station..."
    sleep 10m

done 
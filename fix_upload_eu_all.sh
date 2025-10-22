#!/bin/bash
# Usage: bash fix_upload_eu_all.sh START_DATE END_DATE
# START_DATE and END_DATE in YYYYMMDD format. Script iterates START_DATE .. (END_DATE - 1)
# For each station/date it runs postprocess_eu.py (which handles postprocessing + upload for EU).
# Logs are written to ./logs_eu/

START_DATE=$1
END_DATE=$2

POSTPROCESS_SCRIPT="/flexpart/postprocess_eu.py"
LOG_DIR="./logs_eu"
mkdir -p "${LOG_DIR}"

# Define arguments for each station
ARGS_LIST=(
  "hfd100 0.231 50.977 100 31"
  "htm150 13.4189 56.0976 150 32"
  "gat341 11.4429 53.0657 341 33"
  "kre250 15.08000 49.57200 250 34"
  "lin098 14.12260 52.16630 98 35"
  "ope120 5.50360 48.56190 120 36"
  "ste252 8.45880 53.04310 252 37"
  "trn180 2.11250 47.96470 180 38"
  "hpb131 11.02460 47.80110 131 39"
  "kit200 8.42490 49.09150 200 40"
  "cbw200 4.926161 51.970258 200 41"
  "bik300 23.026750 53.231528 300 42"
  "jfj960 7.9851 46.5475 960 43"
  "oxk163 11.80830 50.03000 163 44"
  "ssl120 7.9113 47.9106 120 45"
  "pal012 24.11570 67.97330 12 46"
  "sac100 2.14200 48.72270 100 47"
  "lmp008 12.63220 35.51810 8 48"
  "pot100 15.7237 40.601 100 49"
  "nor100 17.47940 60.08640 100 50"
  "svb150 19.77500 64.25600 150 51"
)

# basic date validation
date -d "${START_DATE:0:4}-${START_DATE:4:2}-${START_DATE:6:2}" >/dev/null 2>&1 || { echo "Invalid START_DATE"; exit 1; }
date -d "${END_DATE:0:4}-${END_DATE:4:2}-${END_DATE:6:2}" >/dev/null 2>&1 || { echo "Invalid END_DATE"; exit 1; }

current_date=$START_DATE

while [ "$current_date" != "$END_DATE" ]; do
    yyyy=${current_date:0:4}
    mm=${current_date:4:2}
    dd=${current_date:6:2}
    echo "Processing EU date: ${yyyy}-${mm}-${dd}"
    for entry in "${ARGS_LIST[@]}"; do
        read -r station lon lat z cpuid <<< "$entry"
        log_file="${LOG_DIR}/postprocess_eu_${station}_${current_date}.log"
        echo "  Station: ${station}  Date: ${current_date}  CPU: ${cpuid} -> log: ${log_file}"
        # run postprocess.py (assumes it handles EU upload). run in background pinned to cpuid
        nohup bash -c "taskset -c ${cpuid} python3 ${POSTPROCESS_SCRIPT} ${station} ${current_date}" > "${log_file}" 2>&1 &
        sleep 3m
    done

    # increment date by one day
    current_date=$(date -d "${yyyy}-${mm}-${dd} +1 day" +%Y%m%d)
done

echo "All EU jobs launched. Check ${LOG_DIR} for logs."
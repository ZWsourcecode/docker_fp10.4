#!/bin/bash
# Usage: bash fix_upload_global_all.sh START_DATE END_DATE
# START_DATE and END_DATE in YYYYMMDD format. Script will iterate START_DATE .. (END_DATE - 1)
# For each station/date it runs postprocess.py (which handles postprocessing + upload).
# Logs are written to ./logs/

START_DATE=$1
END_DATE=$2

POSTPROCESS_SCRIPT="/flexpart/postprocess.py"
LOG_DIR="./logs"
mkdir -p "${LOG_DIR}"

# Species/station list: "station lon lat z cpuid"
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

# validate dates (basic)
date -d "${START_DATE:0:4}-${START_DATE:4:2}-${START_DATE:6:2}" >/dev/null 2>&1 || { echo "Invalid START_DATE"; exit 1; }
date -d "${END_DATE:0:4}-${END_DATE:4:2}-${END_DATE:6:2}" >/dev/null 2>&1 || { echo "Invalid END_DATE"; exit 1; }

current_date=$START_DATE

while [ "$current_date" != "$END_DATE" ]; do
    yyyy=${current_date:0:4}
    mm=${current_date:4:2}
    dd=${current_date:6:2}
    echo "Processing date: ${yyyy}-${mm}-${dd}"
    for entry in "${ARGS_LIST[@]}"; do
        read -r station lon lat z cpuid <<< "$entry"
        log_file="${LOG_DIR}/postprocess_${station}_${current_date}.log"
        echo "  Station: ${station}  Date: ${current_date}  CPU: ${cpuid} -> log: ${log_file}"
        # run postprocess.py (assumes it handles upload). run in background pinned to cpuid
        nohup bash -c "taskset -c ${cpuid} python3 ${POSTPROCESS_SCRIPT} ${station} ${current_date}" > "${log_file}" 2>&1 &
        # small delay to avoid overwhelming system
        sleep 3m
    done

    # increment date by one day
    current_date=$(date -d "${yyyy}-${mm}-${dd} +1 day" +%Y%m%d)
done

echo "All jobs launched. Check ${LOG_DIR} for logs."
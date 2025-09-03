#!/bin/bash
### run : bash check_missing_csv.sh station_id yyyymm domain
### zhendong.wu@nateko.lu.se

station_id=$1
yyyymm=$2
domain=$3

year=${yyyymm:0:4}
month=${yyyymm:4:2}

if [ "$domain" = "eu" ]; then
    INPUT_DIR="/flexpart/postprocess/flexpartweb/stations_eu/"
else
    INPUT_DIR="/flexpart/postprocess/flexpartweb/stations/"
fi

# Get number of days in the month (portable, no cal)
if [ "${year}${month}" = "$(date +%Y%m)" ]; then
    # If current month, only check up to today
    days_in_month=$(date +%d)
else
    days_in_month=$(date -d "${year}-${month}-01 +1 month -1 day" "+%d")
fi

# Define arguments for each station
station_list=(
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

# station_id_upper=$(echo "${station_id}" | tr '[:lower:]' '[:upper:]')
# TARGET_DIR="${INPUT_DIR}${station_id_upper}/${year}/${month}"
# for day in $(seq -w 1 $days_in_month); do
#     file_name="delta_14C_${station_id}_${yyyymm}${day}.csv"
#     if [ ! -f "${TARGET_DIR}/${file_name}" ]; then
#         echo "Missing: ${file_name}"
#     fi
# done

check_station() {
    local sid="$1"
    local sid_upper
    sid_upper=$(echo "${sid}" | tr '[:lower:]' '[:upper:]')
    local target_dir="${INPUT_DIR}${sid_upper}/${year}/${month}"
    local missing=0
    for day in $(seq -w 1 $days_in_month); do
        if [ "$domain" = "eu" ]; then
            local file_name="delta_14C_${sid}_${yyyymm}${day}eu.csv"
        else
            local file_name="delta_14C_${sid}_${yyyymm}${day}.csv"
        fi
        if [ ! -f "${target_dir}/${file_name}" ]; then
            echo "Missing: ${file_name}"
            missing=1
        fi
    done
    if [ $missing -eq 0 ]; then
        echo "All files are produced for station ${sid_upper} in ${year}${month}"
    fi
}

if [ "$station_id" = "all" ]; then
    for entry in "${station_list[@]}"; do
        sid=$(echo "$entry" | awk '{print $1}')
        check_station "$sid"
    done
else
    check_station "$station_id"
fi

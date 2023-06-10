#!/bin/bash

# location of binary files
BIN="bin"
# location of .gnu files
PLOT="plot"
# name of the executable .gnu file
filename_err="errors"
filename_orb="orbit"

# list of all satellite names available
satellites=("STARLINK" "NUTSAT" "TDRS-3" "ISS")
default_sat="ISS"


make $BIN/main
if [ -z "$1" ] || [[ ! "$satellites" =~ "$1" ]]; then
  echo "No arguments or bad arguments provided. Use the syntax: ./execute.sh <satellite_name>"
  sat_string="Available satellites: "
  for i in "${satellites[@]}"; do
    sat_string+="$i "
  done
  echo $sat_string
  echo "Executing the default command: ./execute.sh $default_sat"
  python3 src/getRV.py $default_sat
  ./$BIN/main $default_sat
  datafile_err="data/${filename_err}/${default_sat}_${filename_err}.txt"
  datafile_orb="data/${filename_orb}/${default_sat}_${filename_orb}.txt"
  datafile_real="data/teme/${default_sat}_1min_interval.txt"
  # datafile="data/sgp4/${default_sat}.txt"
else
  echo "Executing the command: ./execute.sh $1"
  python3 src/getRV.py $1
  ./$BIN/main $1
  datafile_err="data/${filename_err}/${1}_${filename_err}.txt"
  datafile_orb="data/${filename_orb}/${1}_${filename_orb}.txt"
  datafile_real="data/teme/${1}_1min_interval.txt"
  # datafile="data/sgp4/${1}.txt"
fi
echo "Plotting..."
gnuplot -p -c $PLOT/$filename_err.gnu $datafile_err
gnuplot -p -c $PLOT/$filename_orb.gnu $datafile_real $datafile_orb
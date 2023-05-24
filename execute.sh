#!/bin/bash

# location of binary files
BIN="bin"
# location of .gnu files
PLOT="plot"
# name of the executable .gnu file
filename="errors"

# list of all satellite names available
satellites=("STARLINK" "NUTSAT" "TDR-3")


make $BIN/main
if [ -z "$1" ]; then
  echo "No arguments provided. Use the syntax: ./execute.sh <satellite_name>"
  sat_string="Available satellites: "
  for i in "${satellites[@]}"; do
    sat_string+="$i "
  done
  echo $sat_string
  echo "Executing the default command: ./execute.sh STARLINK"
  ./$BIN/main STARLINK
  datafile="data/tle/${filename}_STARLINK.txt"
  echo "Plotting..."
  gnuplot -p -c $PLOT/$filename.gnu $datafile
else
  make $BIN/main
  echo "Executing the command: ./execute.sh $1"
  ./$BIN/main $1
  datafile="data/tle/${filename}_$1.txt"
  echo "Plotting..."
  gnuplot -p -c $PLOT/$filename.gnu $datafile
fi
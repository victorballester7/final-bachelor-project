#!/bin/bash

# location of binary files
BIN="bin"
# location of .gnu files
PLOT="plot"
# name of the executable .gnu file
filename_err="errors"
filename_orb="orbit"

# list of all satellite names available
satellites=("STARLINK_1007" "NUTSAT" "TDRS-3" "ISS" "NAVSTAR_61" "GALILEO_20" "SIRIUS_3" "TESS" "HUBBLE")
true_false=("t" "f")
default_sat="ISS"

name_pointEarth="_pointEarth"
name_sphHarm="_sphHarm"
name_sun="_sun"
name_moon="_moon"
name_solarRad="_solarRad"
name_atmoDrag="_atmoDrag"

make $BIN/main 

if [[ -z "$1" ]] || [[ ! "${satellites[*]}" =~ "$1" ]] || [[ ! "${true_false[*]}" =~ "$2" ]] || [[ ! "${true_false[*]}" =~ "$3" ]] || [[ ! "${true_false[*]}" =~ "$4" ]] || [[ ! "${true_false[*]}" =~ "$5" ]] || [[ ! "${true_false[*]}" =~ "$6" ]] ; then
  echo "No arguments or bad arguments provided. Use the syntax: ./execute.sh <satellite_name> <t/f pointEarth> <t/f sun> <t/f moon> <t/f solarRad> <t/f atmoDrag>"
  sat_string="Available satellites: "
  for i in "${satellites[@]}"; do
    sat_string+="$i "
  done
  echo $sat_string
  echo "Executing the default command: ./execute.sh $default_sat f f f f f"
  satellite=$default_sat
  pointEarth="f"
  sun="f"
  moon="f"
  solarRad="f"
  atmoDrag="f"
else
  echo "Executing the command: ./execute.sh $1"
  satellite=$1
  pointEarth=$2
  sun=$3
  moon=$4
  solarRad=$5
  atmoDrag=$6
fi

end_name=""

if [ "$pointEarth" == "t" ]; then
  end_name+=$name_pointEarth
else
  end_name+=$name_sphHarm
fi
if [ "$sun" == "t" ]; then
  end_name+=$name_sun
fi
if [ "$moon" == "t" ]; then
  end_name+=$name_moon
fi
if [ "$solarRad" == "t" ]; then
  end_name+=$name_solarRad
fi
if [ "$atmoDrag" == "t" ]; then
  end_name+=$name_atmoDrag
fi

end_name+=".txt"

# data files
datafile_err="data/${filename_err}/${satellite}${end_name}"
datafile_orb="data/${filename_orb}/${satellite}${end_name}"
datafile_real="data/teme/${satellite}_1min_interval.txt"

echo "Executing with settings: "
echo "satellite: $satellite"
if [[ "$pointEarth" == "t" ]]; then
  echo "pointEarth: true"
else
  echo "pointEarth: false"
fi
if [[ "$sun" == "t" ]]; then
  echo "sun: true"
else
  echo "sun: false"
fi
if [[ "$moon" == "t" ]]; then
  echo "moon: true"
else
  echo "moon: false"
fi
if [[ "$solarRad" == "t" ]]; then
  echo "solarRad: true"
else
  echo "solarRad: false"
fi
if [[ "$atmoDrag" == "t" ]]; then
  echo "atmoDrag: true"
else
  echo "atmoDrag: false"
fi


echo "Getting the RV data from SGP4..."
python3 src/getRV.py $satellite
echo "Running the integrator..."
./$BIN/main $satellite $pointEarth $sun $moon $solarRad $atmoDrag
echo "Plotting..."
gnuplot -p -c $PLOT/$filename_err.gnu $satellite $datafile_err
# gnuplot -p -c $PLOT/$filename_orb.gnu $datafile_real $datafile_orb
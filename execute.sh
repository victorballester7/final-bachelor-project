#!/bin/bash

# location of binary files
BIN="bin"
# location of .gnu files
PLOT="plot"

if [ -z "$1" ]; then
  echo "No arguments provided. Use the syntax: ./execute.sh <satellite_name>"
  echo "Executing the default command: ./execute.sh ISS"
  ./$BIN/main ISS
  gnuplot -p $PLOT/$ex11.gnu
else
  make $BIN/main
  ./$BIN/main $1
fi

# case $1 in
#   $ex11)
#     make $BIN/$ex11
#     if [ -z "$2" ] && [ -z "$3" ] && [ -z "$4" ]; then
#       echo "No arguments provided. Use the syntax: ./execute.sh $ex11 <hmin> <hmax> <tol>"
#       echo "Executing the default command: ./$BIN/$ex11 0.01 0.05 0.000001"
#       ./$BIN/$ex11 0.01 0.05 0.000001
#     else
#       ./$BIN/$ex11 $2 $3 $4
#     fi
#     gnuplot -p $PLOT/$ex11.gnu
#     ;;
#   $ex14)
#     make $BIN/$ex14
#     echo "Executing $ex14"
#     ./$BIN/$ex14
#     ;;
#   $ex15)
#     make $BIN/$ex15
#     if [ -z "$2" ] && [ -z "$3" ] && [ -z "$4" ] && [ -z "$5" ] && [ -z "$6" ]; then
#       echo "No arguments provided. Use the syntax: ./execute.sh $ex15 <sigma> <rho> <beta> <tf> <nt>"
#       echo "Executing the default command: ./$BIN/$ex15 10 28 2.6 50 5000"
#       ./$BIN/$ex15 10 28 2.6 50 5000
#     else
#       ./$BIN/$ex15 $2 $3 $4 $5 $6
#     fi
#     gnuplot -p $PLOT/$ex15.gnu
#     ;;
#   $ex16)
#     make $BIN/$ex16
#     if [ -z "$2" ]; then
#       echo "No arguments provided. Use the syntax: ./execute.sh $ex16 <mu>"
#       echo "Executing the default command: ./$BIN/$ex16 1.215058560962404e-2"
#       ./$BIN/$ex16 1.215058560962404e-2
#     else
#       ./$BIN/$ex16 $2
#     fi
#     gnuplot -p $PLOT/$ex16.gnu
#     ;;
#   $ex2)
#     make $BIN/$ex2
#     echo "Executing $ex2"
#     ./$BIN/$ex2
#     ;;
#   $ex43)
#     make $BIN/$ex43
#     echo "Executing $ex43"
#     ./$BIN/$ex43
#     ;;
#   $ex44)
#     make $BIN/$ex44
#     if [ -z "$2" ] && [ -z "$3" ] && [ -z "$4" ]; then
#       echo "No arguments provided. Use the syntax: ./execute.sh $ex44 <mu> <tol_newton> <maxit_newton>"
#       echo "Executing the default command: ./$BIN/$ex44 1.215058560962404e-2 1e-12 100"
#       ./$BIN/$ex44 1.215058560962404e-2 1e-12 100
#     else
#       ./$BIN/$ex44 $2 $3 $4
#     fi
#     ;;
#   *)
#     echo "Invalid argument. Use one of the following: $ex11, $ex14, $ex15, $ex16, $ex2, $ex43, $ex44"
#     ;;
# esac


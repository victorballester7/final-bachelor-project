# Get the position and velocity of each TLE from a file
import sys
import os
import numpy as np

if len(sys.argv) != 3:
  print("Usage: python3 getRV.py <file_to_reduce> <times_to_be_reduced>")
  sys.exit(1)

filename = sys.argv[1]
script_dir = os.path.dirname(os.path.abspath(__file__))
filename_in = os.path.join(script_dir + "/../data/errors/" + filename)
filename_out = os.path.join(script_dir + "/../data/plot/" + filename)

with open(filename_in, "r") as f:
  lines = f.readlines()

# print the number of TLEs in filename_out_5
with open(filename_out, "w") as f:
  f.write(lines[0])
  for i in range(1, len(lines), int(sys.argv[2])):
    f.write(lines[i])

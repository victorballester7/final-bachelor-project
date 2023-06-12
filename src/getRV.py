# Get the position and velocity of each TLE from a file
from sgp4.api import Satrec
import sys
import os


if len(sys.argv) != 2:
  print("Usage: python3 getRV.py <satellite_name>")
  sys.exit(1)

satellite_name = sys.argv[1]
script_dir = os.path.dirname(os.path.abspath(__file__))
filename_in = os.path.join(
    script_dir + "/../data/tle/",
    satellite_name + ".txt")
filename_out = os.path.join(
    script_dir + "/../data/teme/",
    satellite_name + "_tles.txt")
filename_out_2 = os.path.join(
    script_dir + "/../data/sgp4/",
    satellite_name + ".txt")
filename_out_3 = os.path.join(
    script_dir + "/../data/teme/",
    satellite_name + "_1h_interval.txt")
filename_out_4 = os.path.join(
    script_dir + "/../data/teme/",
    satellite_name + "_1min_interval.txt")

with open(filename_in, "r") as f:
  lines = f.readlines()

L = []

for i in range(0, len(lines)):
  if i % 2 == 0:
    line1 = lines[i]
    line2 = lines[i + 1]
    satellite = Satrec.twoline2rv(line1, line2)
    jd = satellite.jdsatepoch
    jd_frac = satellite.jdsatepochF
    e, r, v = satellite.sgp4(jd, jd_frac)
    mdjd = jd - 2400000.5
    L.append([mdjd, jd_frac, r, v])

with open(filename_out, "w") as f:
  for i in L:
    f.write(
        str(i[0] + i[1]) + " " +
        str(i[2][0] * 1000) + " " +
        str(i[2][1] * 1000) + " " +
        str(i[2][2] * 1000) + " " +
        str(i[3][0] * 1000) + " " +
        str(i[3][1] * 1000) + " " +
        str(i[3][2] * 1000) + "\n")

# propagate with sgp4 and compare with the original tles
satellite = Satrec.twoline2rv(lines[0], lines[1])
jd = satellite.jdsatepoch
jd_frac = satellite.jdsatepochF
ErrTLE = [[0, 0]]
for i in range(2, len(lines)):
  if i % 2 == 0:
    line1 = lines[i]
    line2 = lines[i + 1]
    satellite_aux = Satrec.twoline2rv(line1, line2)
    T = (satellite_aux.jdsatepoch - jd) + (satellite_aux.jdsatepochF - jd_frac)
    # print(T)
    e, r, v = satellite.sgp4(jd, jd_frac + T)
    e, r0, v0 = satellite_aux.sgp4(
        satellite_aux.jdsatepoch, satellite_aux.jdsatepochF)
    err = sum([(r[i] - r0[i]) ** 2 for i in range(3)])**0.5
    ErrTLE.append([T, err])

with open(filename_out_2, "w") as f:
  for i in ErrTLE:
    f.write(str(i[0]) + " " + str(i[1]) + "\n")


# propagate with sgp4 in a 1h interval
satellite = Satrec.twoline2rv(lines[0], lines[1])
jd = satellite.jdsatepoch
jd_frac = satellite.jdsatepochF
I = []
for i in range(0, 1000):
  e, r, v = satellite.sgp4(jd, jd_frac + i / 24)
  # we do a linear interpolation to get the error from ErrTLE
  elapsed_time = i / (24 * 60)
  err = 0
  for j in range(0, len(ErrTLE) - 1):
    if elapsed_time >= ErrTLE[j][0] and elapsed_time < ErrTLE[j + 1][0]:
      err = (elapsed_time - ErrTLE[j][0]) * (ErrTLE[j + 1][1] -
                                             ErrTLE[j][1]) / (ErrTLE[j + 1][0] - ErrTLE[j][0]) + ErrTLE[j][1]
      break
  I.append([i / 24, r[0], r[1], r[2], v[0], v[1], v[2], e])

mjd = jd - 2400000.5
with open(filename_out_3, "w") as f:
  for i in I:
    f.write(str(mjd + jd_frac + i[0]) + " " +
            str(i[1] * 1000) + " " +
            str(i[2] * 1000) + " " +
            str(i[3] * 1000) + " " +
            str(i[4] * 1000) + " " +
            str(i[5] * 1000) + " " +
            str(i[6] * 1000) + "\n")


# propagate with sgp4 in a 1min interval
satellite = Satrec.twoline2rv(lines[0], lines[1])
jd = satellite.jdsatepoch
jd_frac = satellite.jdsatepochF
I = []
for i in range(0, 20000):
  e, r, v = satellite.sgp4(jd, jd_frac + i / (24 * 60))
  # we do a linear interpolation to get the error from ErrTLE
  elapsed_time = i / (24 * 60)
  err = 0
  for j in range(0, len(ErrTLE) - 1):
    if elapsed_time >= ErrTLE[j][0] and elapsed_time < ErrTLE[j + 1][0]:
      err = (elapsed_time - ErrTLE[j][0]) * (ErrTLE[j + 1][1] -
                                             ErrTLE[j][1]) / (ErrTLE[j + 1][0] - ErrTLE[j][0]) + ErrTLE[j][1]
      break
  I.append([i / (24 * 60), r[0], r[1], r[2], v[0], v[1], v[2], err])

mjd = jd - 2400000.5
with open(filename_out_4, "w") as f:
  for i in I:
    f.write(str(mjd + jd_frac + i[0]) + " " +
            str(i[1] * 1000) + " " +
            str(i[2] * 1000) + " " +
            str(i[3] * 1000) + " " +
            str(i[4] * 1000) + " " +
            str(i[5] * 1000) + " " +
            str(i[6] * 1000) + " " +
            str(i[7]) + "\n")

# print("Done!")

# bstar = satellite.bstar

# print("B* = " + str(bstar))

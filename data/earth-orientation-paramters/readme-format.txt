# from https://maia.usno.navy.mil/products/daily

The format of the finals.data, finals.daily, and finals.all files is:

Col.#    Format  Quantity
-------  ------  -------------------------------------------------------------
1-2      I2      year (to get true calendar year, add 1900 for MJD<=51543 or add 2000 for MJD>=51544)
3-4      I2      month number
5-6      I2      day of month
7        X       [blank]
8-15     F8.2    fractional Modified Julian Date (MJD UTC)
16       X       [blank]
17       A1      IERS (I) or Prediction (P) flag for Bull. A polar motion values
18       X       [blank]
19-27    F9.6    Bull. A PM-x (sec. of arc)
28-36    F9.6    error in PM-x (sec. of arc)
37       X       [blank]
38-46    F9.6    Bull. A PM-y (sec. of arc)
47-55    F9.6    error in PM-y (sec. of arc)
56-57    2X      [blanks]
58       A1      IERS (I) or Prediction (P) flag for Bull. A UT1-UTC values
59-68    F10.7   Bull. A UT1-UTC (sec. of time)
69-78    F10.7   error in UT1-UTC (sec. of time)
79       X       [blank]
80-86    F7.4    Bull. A LOD (msec. of time) -- NOT ALWAYS FILLED
87-93    F7.4    error in LOD (msec. of time) -- NOT ALWAYS FILLED
94-95    2X      [blanks]
96       A1      IERS (I) or Prediction (P) flag for Bull. A nutation values
97       X       [blank]
98-106   F9.3    Bull. A dPSI (msec. of arc)
107-115  F9.3    error in dPSI (msec. of arc)
116      X       [blank]
117-125  F9.3    Bull. A dEPSILON (msec. of arc)
126-134  F9.3    error in dEPSILON (msec. of arc)
135-144  F10.6   Bull. B PM-x (sec. of arc)
145-154  F10.6   Bull. B PM-y (sec. of arc)
155-165  F11.7   Bull. B UT1-UTC (sec. of time)
166-175  F10.3   Bull. B dPSI (msec. of arc)
176-185  F10.3   Bull. B dEPSILON (msec. of arc)


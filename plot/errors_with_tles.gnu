# Set the output file format and name
set encoding utf8; 

# Set the title and labels for the graph
set title 'Error between TLE and integrated orbit of the ' . ARG1 . ' satellite'
set xlabel 'Days'
set ylabel 'Position error'


# define the data file and block separator
datafile = ARG2
datafile2 = ARG3

stats datafile using 1 nooutput
maxDay=ceil(STATS_max) # round up to nearest integer of the maximum temperature in the data file (which is in the 3rd column)

# Set the range for the X and Y axes
set xrange [*:maxDay]
set yrange [*:*]


# Set the style of the data points and lines
set style data points
set style line 1 lc rgb '#ff0000' pt 7 ps 0.5 lw 2


# Plot the data from the datafile
plot datafile using 1:2 smooth mcs with lines, datafile2 using 1:2 smooth mcs with lines
# plot datafile using 1:2:3 with yerrorbars ls 1 lc 1
# plot datafile using 1:2 with linespoints ls 1 lc 1

# Uncomment the line below if you want to display the graph window
# pause -1

# Set the output file format and name
set encoding utf8; 

# Set the title and labels for the graph
set title 'Error between TLE and integrated orbit of the ' . ARG1 . ' satellite'
set xlabel 'Days'
set ylabel 'Position error'

# Set the range for the X and Y axes
set xrange [*:*]
set yrange [*:*]

# define the data file and block separator
datafile = ARG2

# Set the style of the data points and lines
set style data points
set style line 1 lc rgb '#0060ad' pt 7 ps 0.5 lw 2
set style fill transparent solid 0.15 noborder


# Plot the data from the datafile
plot datafile using 1:($2-$3):($2+$3) with filledcurves, datafile using 1:2 smooth mcs with lines
# plot datafile using 1:2:3 with yerrorbars ls 1 lc 1
# plot datafile using 1:2 with linespoints ls 1 lc 1

# Uncomment the line below if you want to display the graph window
# pause -1

# Set the output file format and name
set encoding utf8; 

# Set the title and labels for the graph
set title '2D Graph'
set xlabel 'Days'
set ylabel 'Position error'

# Set the range for the X and Y axes
set xrange [*:*]
set yrange [*:*]

# define the data file and block separator
datafile = ARG1
# datafile = 'data/tle/errors_STARLINK.txt'

# Set the style of the data points and lines
set style data points
set style line 1 lc rgb '#0060ad' pt 7 ps 0.5

# Plot the data from the datafile
plot datafile using 1:2 with linespoints ls 1 lc 1 title "Error in L1"

# Uncomment the line below if you want to display the graph window
# pause -1

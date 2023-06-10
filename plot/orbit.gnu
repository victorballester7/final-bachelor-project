set encoding utf8; 
# set terminal qt size 1750, 850;

# set the colors for each block
set linetype 1 linecolor rgb "red"
set linetype 2 linecolor rgb "blue"
set linetype 3 linecolor rgb "green"
set linetype 4 linecolor rgb "orange"
set linetype 5 linecolor rgb "purple"
set linetype 6 linecolor rgb "gray"
set linetype 7 linecolor rgb "black"

set view equal xyz
set xlabel 'x'; 
set ylabel 'y'; 
set zlabel 'z';
set key outside bottom center; 

# define the data file and block separator
datafile_real = ARG1
datafile_integ = ARG2
# datafile = 'data/tle/orbit_STARLINK.txt'


# count the number of blocks in the file
# stats datafile_integ using 2:3 nooutput
# nblocks = STATS_blocks

# plot each block with a for loop
# plotcmd = ''
# do for [i=0:nblocks-1] {
#   title = system(sprintf("sed -n '1{p;q;}' %s | cut -d ' ' -f 2-", datafile))
#   title = sprintf("'%s'", title)
#   plotcmd = sprintf("%s '%s' index %d using 2:3 with lines title %s linestyle %d,", plotcmd, datafile, i, title, i+1)
#   if (i < nblocks-1) { plotcmd = sprintf("%s %s", plotcmd, separator) }
# }
# plot plotcmd
# splot datafile index 0 using 1:2:3 title columnheader(1) pointtype 7 pointsize 1 linecolor 1, \
#       for[i=1:nblocks-1] datafile index i using 1:2:3  with lines linecolor i+1
# splot for[i=1:nblocks-1] datafile index i using 1:2:3  with lines linecolor i+1

splot datafile_real using 2:3:4 with lines linecolor 1,datafile_integ using 1:2:3 with lines linecolor 2
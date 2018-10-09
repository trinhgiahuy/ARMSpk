set terminal svg size 1600,600 dynamic enhanced fname 'Times' fsize 28 butt dashlength 1.0
set output "../figures/sys-util.svg"

set auto x
set style data histogram
#set style histogram rowstacked gap 1
set style histogram rowstacked title textcolor lt -1
set style data histograms
#set style fill solid border -1
set style fill   solid 1.00 border lt -1
set boxwidth 0.8
#set key Left reverse box
set key outside right top vertical Left reverse noenhanced autotitle columnhead nobox
set key invert samplen 4 spacing 1 width 0 height 0
set xtic font ",24" rotate by -45 scale 0 left
set yrange [0:100]
set ylabel "HPC resource utilization [%]"
set grid

plot \
	"../data/sys-util.data" using  2:xtic(1) title 'geo', \
	"" using  3 title 'chm', \
	"" using  4 title 'phy', \
	"" using  5 title 'qcd', \
	"" using  6 title 'mat', \
	"" using  7 title 'eng', \
	"" using  8 title 'mcs', \
	"" using  9 title 'bio', \
	"" using 10 title 'oth' fs pattern 1 fs solid 0.25, \


set terminal svg size 1600,500 dynamic enhanced fname 'Times' fsize 24 butt dashlength 1.0
set output "../figures/alu-fpu-ops.svg"
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
set xtic font ",22" rotate by -45 scale 0 left
set yrange [0:100]
set ylabel "Percentage of Operations [%]"
set grid

plot \
	"../data/alu-fpu-ops.data" using  3:xtic(1) title 'FP64' lc rgb "#AAC864", \
	"" using  5 title 'FP32' lc rgb "#7A1420", \
	"" using  7 title 'INT' lc rgb "#6E92A1"


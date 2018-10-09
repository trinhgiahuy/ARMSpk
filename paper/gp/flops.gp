set terminal svg size 1600,600 dynamic enhanced fname 'Times' fsize 28 butt dashlength 1.0
set output "../figures/flops.svg"

set grid
set auto x
set auto y
set xrange [0:24]
set yrange [0.1:1e4]
set xtic font ",24" rotate by -45 scale 0 left
set logscal y

set key right bottom vertical Right box width +2
#reverse noenhanced autotitle columnhead box
set datafile missing '-'

bdw = "#A61A00"
knl = "#00B358"
knm = "#1924B1"

set ylabel "FP (SP+DP) Operations per Second [Gflop/s]"

# max DP flops lines
bdwmax(x)=( -0.5 < x && x < 24.5 ) ? 345.6 : 1/0
knlmax(x)=( -0.5 < x && x < 24.5 ) ? 2662 : 1/0
knmmax(x)=( -0.5 < x && x < 24.5 ) ? 1728 : 1/0

plot \
	"../data/flops.data" u 2:xtic(1) pt 3 ps 0.8 lc rgb bdw title 'BDW', \
	bdwmax(x) with lines lt 0 lw 2 lc rgb bdw notitle, \
	"" u 4:xtic(1) pt 20  ps 0.8 lc rgb knl title 'KNL', \
	knlmax(x) with lines lt 0 lw 2 lc rgb knl notitle, \
	"" u 6:xtic(1) pt 9  ps 0.8 lc rgb knm title 'KNM', \
	knmmax(x) with lines lt 0 lw 2 lc rgb knm notitle


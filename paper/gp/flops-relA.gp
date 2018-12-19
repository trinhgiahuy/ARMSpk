set terminal svg size 1600,600 dynamic enhanced fname 'Times' fsize 28 butt dashlength 1.0
set output "../figures/flops-relA.svg"

set grid
set auto x
set auto y
set xrange [-0.5:23.5]
set yrange [0:4]
set ytics 0,1,4
set xtic font ",24" rotate by -45 scale 0 left

set key left top vertical Right maxrows 3 box width +2
set datafile missing '-'

bdw = "#A61A00"
knl = "#00B358"
knm = "#1924B1"

set ylabel "Rel. Perf. (Gflop/s) Improvement over BDW"

# max DP flops lines
bdwmax(x)=( -1 < x && x < 24 ) ? 1.0 : 1/0

plot \
	"../data/flops.data" u ($4/$2):xtic(1) pt 20  ps 0.8 lc rgb knl title 'KNL_{rel}', \
	"" u ($6/$2):xtic(1) pt 9  ps 0.8 lc rgb knm title 'KNM_{rel}', \
	bdwmax(x) with lines lt 0 lw 2 lc rgb bdw title 'BDW_{rel}'


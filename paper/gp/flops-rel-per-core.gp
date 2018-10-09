set terminal svg size 1600,600 dynamic enhanced fname 'Times' fsize 28 butt dashlength 1.0
set output "../figures/flops-rel-per-core.svg"

set grid
set auto x
set auto y
set xrange [-0.5:23.5]
set yrange [0:2]
set ytics 0,1,2
set xtic font ",24" rotate by -45 scale 0 left

set key right top vertical Right box width +2
#reverse noenhanced autotitle columnhead box
set datafile missing '-'

bdw = "#A61A00"
knl = "#00B358"
knm = "#1924B1"

set ylabel "Rel. Perf. (Gflop/s) Improv. (per core)"

# max DP flops lines
bdwmax(x)=( -1 < x && x < 24 ) ? 1.0 : 1/0

plot \
	"../data/flops.data" u (($4/$5)/($2/$3)):xtic(1) pt 20  ps 0.8 lc rgb knl title 'KNL', \
	"" u (($6/$7)/($2/$3)):xtic(1) pt 9  ps 0.8 lc rgb knm title 'KNM', \
	bdwmax(x) with lines lt 0 lw 2 lc rgb bdw title 'BDW'


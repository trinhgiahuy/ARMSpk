set terminal svg size 1600,600 dynamic enhanced fname 'Times' fsize 28 butt dashlength 1.0
set output "../figures/t2s-rel.svg"

set grid
set auto x
set auto y
set xrange [-0.5:23.5]
set yrange [0.5:1.5]
set ytics 0.5,0.5,1.5
set xtic font ",24" rotate by -45 scale 0 left

set key left top vertical Right box width +2
#reverse noenhanced autotitle columnhead box
set datafile missing '-'

bdw = "#A61A00"
knl = "#00B358"
knm = "#1924B1"

set ylabel "Speedup (w.r.t Time-to-Solution)"

knlmax(x)=( -1 < x && x < 24 ) ? 1.0 : 1/0

plot \
	"../data/t2solv.data" u ($4/$6):xtic(1) pt 20  ps 0.8 lc rgb knm title 'KNM', \
	knlmax(x) with lines lt 0 lw 2 lc rgb knl title 'KNL'


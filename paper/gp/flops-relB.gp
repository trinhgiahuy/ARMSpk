set terminal svg size 1600,600 dynamic enhanced fname 'Times' fsize 28 butt dashlength 1.0
set output "../figures/flops-relB.svg"

set grid
set auto x
set auto y
set xrange [-0.5:23.5]
set yrange [0:100]
set ytics 0,20,100
set xtic font ",24" rotate by -45 scale 0 left

set key left top vertical Right maxrows 3 box width +2
set datafile missing '-'

bdw = "#A61A00"
knl = "#00B358"
knm = "#1924B1"

set ylabel "Abs. achieved Gflop/s out of Peak [in %]"

# max DP flops lines
bdwmax(x)=( -1 < x && x < 24 ) ? 1.0 : 1/0

plot \
	"../data/flops.data" u (100.0*$4/$9):xtic(1)  pt 20 ps 0.8 lc rgb knl title 'KNL_{abs}', \
	"" u (100.0*$6/$10) pt 9 ps 0.8 lc rgb knm title 'KNM_{abs}', \
	"" u (100.0*$2/$8)  pt 3 ps 0.8 lc rgb bdw title 'BDW_{abs}'

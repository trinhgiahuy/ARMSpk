set terminal svg size 1600,600 dynamic enhanced fname 'Times' fsize 28 butt dashlength 1.0
set output "../figures/flops-rel.svg"

set grid
set auto x
set auto y
set xrange [-0.5:23.5]
set yrange [0:4]
set y2range [0:100]
set ytics 0,1,4
set y2tics 0,20,100
set xtic font ",24" rotate by -45 scale 0 left

set key left top vertical Right maxrows 3 box width +2
set datafile missing '-'

bdw = "#A61A00"
knl = "#00B358"
knm = "#1924B1"
bdw2 = "#D18A7D"
knl2 = "#7DD8A9"
knm2 = "#898FD7"

set ylabel "Rel. Perf. (Gflop/s) Improvement over BDW"
set y2label "Abs. achieved Gflop/s out of Peak [in %]"

# max DP flops lines
bdwmax(x)=( -1 < x && x < 24 ) ? 1.0 : 1/0

plot \
	"../data/flops.data" u ($4/$2):xtic(1) pt 20  ps 0.8 lc rgb knl title 'KNL_{rel}', \
	"" u ($6/$2):xtic(1) pt 9  ps 0.8 lc rgb knm title 'KNM_{rel}', \
	bdwmax(x) with lines lt 0 lw 2 lc rgb bdw title 'BDW_{rel}', \
	"" u (100.0*$4/$9)  pt 7 ps 0.6 lc rgb knl2 title 'KNL_{abs}' axes x1y2 , \
	"" u (100.0*$6/$10) pt 7 ps 0.6 lc rgb knm2 title 'KNM_{abs}' axes x1y2, \
	"" u (100.0*$2/$8)  pt 7 ps 0.6 lc rgb bdw2 title 'BDW_{abs}' axes x1y2

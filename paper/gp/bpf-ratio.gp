set terminal svg size 1600,600 dynamic enhanced fname 'Times' fsize 28 butt dashlength 1.0
set output "../figures/bpf-ratio.svg"

set grid
set auto x
set auto y
set xrange [0.1:24.9]
set yrange [1e-2:700]
set xtic font ",24" rotate by -45 scale 0 left
set logscale y
set format y "10^{%L}"

set key right top vertical Right box width +2
#reverse noenhanced autotitle columnhead box
set datafile missing '-'

bdw = "#A61A00"
knl = "#00B358"
knm = "#1924B1"

set ylabel "Bytes Transfered per FP (SP+DP) Operation"

bdwmax(x)=( -0.5 < x && x < 25.4 ) ? 122.0/691.0 : 1/0
knlmax(x)=( -0.5 < x && x < 25.4 ) ? 439.0/5324.0 : 1/0
knmmax(x)=( -0.5 < x && x < 25.4 ) ? 430.0/13824.0 : 1/0

plot \
	"../data/bytes-n-flops.data" u (($2*$5)/$3):xtic(1) pt 3 ps 0.8 lc rgb bdw title 'BDW', \
	bdwmax(x) with lines lt 0 lw 2 lc rgb bdw notitle, \
	"" u    (($6*$9)/$7) pt 20  ps 0.8 lc rgb knl title 'KNL', \
	knlmax(x) with lines lt 0 lw 2 lc rgb knl notitle, \
	"" u (($10*$13)/$11) pt 9  ps 0.8 lc rgb knm title 'KNM', \
	knmmax(x) with lines lt 0 lw 2 lc rgb knm notitle


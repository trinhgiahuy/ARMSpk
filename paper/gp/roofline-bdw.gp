set terminal svg size 1600,1200 dynamic enhanced fname 'Times' fsize 32 butt dashlength 1.0
set output "../figures/roofline-bdw.svg"

# gflops
knl_fpeak = 2662.0
knm_fpeak = 1728.0
bdw_fpeak = 691.0
# gb/s
knl_mpeak = 439.0
knm_mpeak = 430.0
bdw_mpeak = 122.0

xmin = 0.001
xmax = 100
ymin = 0.1
ymax = 2000
set xtics nomirror
set xrange [xmin:xmax]
set logscale x 10
set yrange [ymin:ymax]
set logscale y 10

#	Functions
mem(x,y)     = exp( log( y ) - log( x ))
min(a,b) = (a < b) ? a : b
max(a,b) = (a > b) ? a : b
knl_froof(x) = knl_fpeak
knl_mroof(x) = mem(knl_fpeak / knl_mpeak, knl_fpeak) * x
knl_rflne(x) = min(knl_froof(x), knl_mroof(x))
knm_froof(x) = knm_fpeak
knm_mroof(x) = mem(knm_fpeak / knm_mpeak, knm_fpeak) * x
knm_rflne(x) = min(knm_froof(x), knm_mroof(x))
bdw_froof(x) = bdw_fpeak
bdw_mroof(x) = mem(bdw_fpeak / bdw_mpeak, bdw_fpeak) * x
bdw_rflne(x) = min(bdw_froof(x), bdw_mroof(x))

set grid
#set key left top vertical Right box width +2
unset key
set xlabel "Arithmetic Intensity (flop/byte)"
set ylabel "Gflop/s"

bdw = "#A61A00"
knl = "#00B358"
knm = "#1924B1"

set label 1 "Theor. Peak Performance (FP64)" at xmax-10, 1.25*bdw_froof(xmax) right
set label 2 "Stream Triad Bandwidth (GB/s)" at 1.25*xmin, 1.6*bdw_mroof(xmin) left rotate by 42

plot bdw_rflne(x) lt 1 lc rgb "black" lw 4 notitle, \
     "../data/bytes-n-flops.data" u ($3/$5)/($2):($3/$5) pt 28 ps 0.6 lc rgb bdw title 'BDW', \
     "" u ($3/$5)/($2):($3/$5):($1) with labels offset 0,-1 font "Times,22" point pt 28 ps 0.6 lc rgb bdw notitle

#plot knl_rflne(x) lt 1 lc rgb knl lw 4 notitle, \
#     knm_rflne(x) lt 1 lc rgb knm lw 4 notitle, \
#     bdw_rflne(x) lt 1 lc rgb bdw lw 4 notitle, \
#     "../data/bytes-n-flops.data" u ($7/$9)/($6):($7/$9) pt 20 ps 0.6 lc rgb knl title 'KNL', \
#     "" u ($11/$13)/($10):($11/$13) pt 9 ps 0.6 lc rgb knm title 'KNM', \
#     "" u ($11/$13)/($10):($11/$13):($1) with labels offset -2.5,-.3 font "Times,24" point pt 9 ps 0.6 lc rgb knm notitle, \
#     "" u ($3/$5)/($2):($3/$5) pt 28 ps 0.6 lc rgb bdw title 'BDW', \
#     "" u ($3/$5)/($2):($3/$5):($1) with labels offset 0,-1 font "Times,24" point pt 28 ps 0.6 lc rgb bdw notitle


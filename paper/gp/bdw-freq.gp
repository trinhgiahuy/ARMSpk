set terminal svg size 1600,500 dynamic enhanced fname 'Times' fsize 22 butt dashlength 1.0
set output "../figures/bdw-freq.svg"
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
set xtic font ",20" rotate by -45 scale 0 left
set yrange [1:2]
set ylabel "Speedup"
set grid
set label "Broadwell-EP (2x)" left at -0.5,1.95

plot \
	"../data/bdw-freq.data" using ($2<0 ? NaN : $2):xtic(1) title '1.2 GHz', \
	"" using ($3<0 ? NaN : $3) title '1.3 GHz', \
	"" using ($4<0 ? NaN : $4) title '1.4 GHz', \
	"" using ($5<0 ? NaN : $5) title '1.5 GHz', \
	"" using ($6<0 ? NaN : $6) title '1.6 GHz', \
	"" using ($7<0 ? NaN : $7) title '1.7 GHz', \
	"" using ($8<0 ? NaN : $8) title '1.8 GHz', \
	"" using ($9<0 ? NaN : $9) title '1.9 GHz' fs pattern 1 fs solid 0.5, \
	"" using ($10<0 ? NaN : $10) title '2.0 GHz' fs pattern 1 fs solid 0.5, \
	"" using ($11<0 ? NaN : $11) title '2.1 GHz' fs pattern 1 fs solid 0.5, \
	"" using ($12<0 ? NaN : $12) title '2.2 GHz' fs pattern 1 fs solid 0.5, \
	"" using ($13<0 ? NaN : $13) title '2.2 GHz +TB' fs pattern 1 fs solid 0.25


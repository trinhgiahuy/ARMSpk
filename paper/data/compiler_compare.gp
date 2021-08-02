set term pdf color enhanced font "arial,11" size 7.5,16
set output "compiler_compare.pdf"

unset key
set style increment default
set view map scale 1
set style data lines
set xtics border in scale 0,0 mirror norotate  autojustify
set ytics border in scale 0,0 mirror norotate  autojustify
set ztics border in scale 0,0 nomirror norotate  autojustify
unset cbtics
set rtics axis in scale 0,0 nomirror norotate  autojustify
set xlabel "Compiler Variant" font "arial,13"
set xrange [-0.5:4.5] noreverse nowriteback
set yrange [-0.5:115.5] noreverse nowriteback
set x2range [*:*] noreverse writeback
set y2range [*:*] noreverse writeback
set zrange [*:*] noreverse writeback

set cblabel "Relative Performance Gain" font "arial,13"
set cbrange [-1:1] noreverse nowriteback
set rrange [*:*] noreverse writeback
set palette defined ( 0 "#DE77AE", 1 "#E9A3C9", 2 "#ffffff", 3 "#7FBC41", 4 "#4DAC26" )

set datafile separator ","

set rmargin 0

set colorbox

set label front "Micro Kernels (all with [1|12])" at 2.2,115.0 font "arial,13"
set arrow from -0.5,114.5 to 4.5,114.5 nohead lc rgb "#000000" front

set arrow from -0.5,92.5 to 4.5,92.5 nohead lc rgb "#000000" front
set label front "PolyBench (all with [1|1])" at 2.2,92.0 font "arial,13"
set arrow from -0.5,91.5 to 4.5,91.5 nohead lc rgb "#000000" front

set arrow from -0.5,61.5 to 4.5,61.5 nohead lc rgb "#000000" front
set label front "Ranking" at 2.2,61.0 font "arial,13"
set arrow from -0.5,60.5 to 4.5,60.5 nohead lc rgb "#000000" front

set arrow from -0.5,56.5 to 4.5,56.5 nohead lc rgb "#000000" front
set label front "ECP proxy apps" at 2.2,56.0 font "arial,13"
set arrow from -0.5,55.5 to 4.5,55.5 nohead lc rgb "#000000" front

set arrow from -0.5,44.5 to 4.5,44.5 nohead lc rgb "#000000" front
set label front "RIKEN miniapps" at 2.2,44.0 font "arial,13"
set arrow from -0.5,43.5 to 4.5,43.5 nohead lc rgb "#000000" front

set arrow from -0.5,36.5 to 4.5,36.5 nohead lc rgb "#000000" front
set label front "SPEC CPU int (all with [1|1])" at 2.2,36.0 font "arial,13"
set arrow from -0.5,35.5 to 4.5,35.5 nohead lc rgb "#000000" front

set arrow from -0.5,25.5 to 4.5,25.5 nohead lc rgb "#000000" front
set label front "SPEC CPU float" at 2.2,25.0 font "arial,13"
set arrow from -0.5,24.5 to 4.5,24.5 nohead lc rgb "#000000" front

set arrow from -0.5,14.5 to 4.5,14.5 nohead lc rgb "#000000" front
set label front "SPEC OMP" at 2.2,14.0 font "arial,13"
set arrow from -0.5,13.5 to 4.5,13.5 nohead lc rgb "#000000" front

set label front "-1.0" at 4.7,-1.1 font "arial,13"
set label front "+1.0" at 4.7,116.1 font "arial,13"
set title "Time-to-Solution [in s] (FJtrad) and Relative Performance Gain (others)" font "arial,16" offset 0,-0.5

set arrow from .5,-.5 to .5,115.5 nohead lc rgb "#000000" front
set arrow from .48,-.5 to .48,115.5 nohead lc rgb "#000000" front
set label front "runtime" at -0.2,115.0 font "arial,13"

plot "compiler_compare.data" \
        using (int(int($1)/116)):(int($1)%116):5:xtic(3):ytic(2) with image, \
     "" using (int(int($1)/116)):(int($1)%116):7 with labels
#plot "compiler_compare.data" \
#        using (int(int($1)/116)):(int($1)%116):5:xtic(3):ytic(2) with image, \
#     "" using (int(int($1)/116)):(int($1)%116):(sprintf("%.3f",($4))) with labels

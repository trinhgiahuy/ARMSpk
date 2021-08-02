set term pdf color enhanced font "arial,11" size 4.5,2.6
set output "demo_polly.pdf"

# Don't show the legend in the chart
set nokey

# Thinner, filled bars
set boxwidth 0.5
set style fill solid 1.00

# Set a title and Y label. The X label is obviously months, so we don't set it.
set title "Polyhedral optimizations for Xeon (Intel icc) vs. A64FX (FJtrad)" font ",14" tc rgb "#606060"
set ylabel "Speedup advantage of Intel icc"

# Rotate X labels and get rid of the small stripes at the top (nomirror)
set xtics nomirror rotate by -45

# Show human-readable Y-axis. E.g. "100 k" instead of 100000.
set format y '%.0sx'

# Replace small stripes on the Y-axis with a horizontal gridlines
set tic scale 0
set grid ytics lc rgb "#505050"

# Remove border around chart
unset border

# Manual set the Y-axis range
set yrange [0:]

set datafile separator ","

plot "demo_polly.data" using 1:2:xticlabels(3) with boxes lt rgb "#406090"

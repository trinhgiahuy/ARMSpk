#!/usr/bin/env python

from common import *

figure_setup("Throughput", "GB/sec", (10, 5), 50)

bw_app = []
bw_xeon = []
bw_knl = []
bw_knm = []
for app, rows in RESULT_DIC.items():
    if not is_plotted(rows):
        continue
    bw_app.append(app)
    bw_xeon.append(get_bandwidth(rows[TITLE_XEON]))
    bw_knl.append(get_bandwidth(rows[TITLE_KNL]))
    bw_knm.append(get_bandwidth(rows[TITLE_KNM]))

plt.axhline(y=PEAKBW_XEON, color=COLOR_XEON, ls='--')
plt.plot(bw_app, bw_xeon, ".", color=COLOR_XEON, label="Xeon", markersize=30)

plt.axhline(y=PEAKBW_KNL, color=COLOR_KNL, ls='--')
plt.plot(bw_app, bw_knl, ".", color=COLOR_KNL, label="KNL", markersize=30)

plt.axhline(y=PEAKBW_KNM, color=COLOR_KNM, ls='--')
plt.plot(bw_app, bw_knm, ".", color=COLOR_KNM, label="KNM", markersize=30)

figure_save("img/throughput.eps")

#!/usr/bin/env python

from common import *

figure_setup("FLOPS", "GFLOP/s", (10, 5), 50)

flops_app = []
flops_xeon = []
flops_knl = []
flops_knm = []
for app, rows in RESULT_DIC.items():
    if not is_plotted(rows):
        continue
    flops_app.append(app)
    flops_xeon.append(get_flops(rows[TITLE_XEON]))
    flops_knl.append(get_flops(rows[TITLE_KNL]))
    flops_knm.append(get_flops(rows[TITLE_KNM]))

plt.plot(flops_app, flops_xeon, ".", color=COLOR_XEON, label="XEON", markersize=30)
plt.plot(flops_app, flops_knl, ".", color=COLOR_KNL, label="KNL", markersize=30)
plt.plot(flops_app, flops_knm, ".", color=COLOR_KNM, label="KNM", markersize=30)

figure_save("img/flops.eps")

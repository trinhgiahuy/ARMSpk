#!/usr/bin/env python

from common import *

figure_setup("Ratio", "%", (12, 5), 20)

ratio_app = range((len(RESULT_DIC) - 4) * 5)
ratio_label = []
ratio_fp64 = []
ratio_fp32 = []
ratio_int = []
for app, rows in RESULT_DIC.items():
    if not is_plotted(rows):
        continue

    ratio_label.append("")
    ratio_label.append("")
    ratio_label.append(app)
    ratio_label.append("")
    ratio_label.append("")

    ratio_xeon = get_ratio(rows[TITLE_XEON])
    ratio_knl  = get_ratio(rows[TITLE_KNL])
    ratio_knm  = get_ratio(rows[TITLE_KNM])

    ratio_fp64.append(BOTTOM)
    ratio_fp64.append(ratio_xeon[0])
    ratio_fp64.append(ratio_knl[0])
    ratio_fp64.append(ratio_knm[0])
    ratio_fp64.append(BOTTOM)

    ratio_fp32.append(BOTTOM)
    ratio_fp32.append(ratio_xeon[0] + ratio_xeon[1])
    ratio_fp32.append(ratio_knl[0] + ratio_knl[1])
    ratio_fp32.append(ratio_knm[0] + ratio_knm[1])
    ratio_fp32.append(BOTTOM)

    ratio_int.append(BOTTOM)
    ratio_int.append(100)
    ratio_int.append(100)
    ratio_int.append(100)
    ratio_int.append(BOTTOM)

plt.xticks(np.arange(len(ratio_app)), ratio_label)

plt.bar(ratio_app, ratio_int,  0.8, label="INT")
plt.bar(ratio_app, ratio_fp32, 0.8, label="FP32")
plt.bar(ratio_app, ratio_fp64, 0.8, label="FP64")

figure_save("img/opratio.eps")

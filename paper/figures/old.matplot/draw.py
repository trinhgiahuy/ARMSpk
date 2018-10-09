import matplotlib
matplotlib.use("Agg")

import numpy as np
import warnings
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns

def figure_setup(title, ylabel, figsize, grid_interval):
    plt.figure(figsize=figsize)
    sns.set('paper', 'whitegrid', 'deep',
            font_scale=1.5,
            rc={"lines.linewidth": 2,
                'grid.linestyle': '--'})
    plt.title(title, fontsize=20, fontweight='bold')
    plt.ylabel(ylabel)
    plt.axes().yaxis.set_major_locator(ticker.MultipleLocator(grid_interval))

def figure_save(filename):
    plt.ylim(ymin=0)
    plt.xticks(rotation=290)
    plt.legend(loc='center left',
               bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    plt.savefig(filename)

figure_setup("TEMP", "", (10,10), 1) # To change the color or something
colors=plt.rcParams['axes.prop_cycle'].by_key()['color']
COLOR_XEON=colors[2]
COLOR_KNL=colors[0]
COLOR_KNM=colors[1]

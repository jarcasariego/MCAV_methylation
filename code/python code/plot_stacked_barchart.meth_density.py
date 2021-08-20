#!/usr/bin/env python3

"""
> plot_stacked_barchart.meth_density.py <

Plots stacked barchart to illustrate the changes in frequencies of repeat
elements across temperatures.
"""
import csv
from pathlib import Path

from matplotlib import cm
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats
import seaborn as sns

# read data
df = []

tsv_files = Path.cwd().glob('???????.meth_density.tsv')
for t in tsv_files:
    temp = pd.read_table(t, usecols=[0, 11, 15, 19])
    temp['treatment'] = t.name[:2]
    temp['replicate'] = t.name[2:6]
    df.append(temp)

df = pd.concat(df, axis=0, ignore_index=True)
df.columns = [x.replace('within_meth_', '') for x in df.columns]

# exclude sums and nots, because this plot focuses on
# individual repeat element types
excluded_repeats = ['RC','Satellite', 'Unspecified', 'snRNA', 'sum_repeats', 'not_repeats']
df = df[~df['repeat_family'].isin(excluded_repeats)]
df = df.reset_index(drop=True)

# harness the power of "groupby"s to calculate stats of specific groups
df_means = df.groupby(['repeat_family', 'treatment']).mean()
df_sems = df.drop(['replicate'], axis=1).groupby(['repeat_family', 'treatment']).sem()

# convert fractions to percentages for plotting
df_means = df_means * 100
df_sems = df_sems * 100

# start plotting
repeat_types = df_means.index.get_level_values(0).unique()
yloc = [0.2, 0.5, 0.8, 1.1, 1.5, 1.8, 2.1, 2.4, 2.8, 3.1, 3.4, 3.7]
col = cm.get_cmap('rainbow')

sns.set_style('ticks')
f, ax = plt.subplots(figsize=(10, 5))

barh_plots = {}
for n, r in enumerate(repeat_types):
    barh_col = col((n + 0.5) / len(repeat_types))
    if n == 0:
        barh_plots[r] = plt.barh(yloc, df_means.loc[r].unstack(), 0.3,
                                 color=barh_col, ecolor=barh_col,
                                 xerr=df_sems.loc[r].unstack())
        baseline = df_means.loc[r].unstack()
    else:
        barh_plots[r] = plt.barh(yloc, df_means.loc[r].unstack(), 0.3,
                                 color=barh_col, ecolor=barh_col,
                                 left=baseline,
                                 xerr=df_sems.loc[r].unstack())
        baseline += df_means.loc[r].unstack()

# annotations
plt.xlabel('% methylated positions in repeats')
plt.yticks(yloc, ['BC', 'BH', 'CC', 'CH'] * 3)
plt.text(-0.5, 0.5, 'CpG', verticalalignment='center', rotation=90)
plt.text(-0.5, 1.8, 'CHG', verticalalignment='center', rotation=90)
plt.text(-0.5, 3.1, 'CHH', verticalalignment='center', rotation=90)
plt.legend([barh_plots[x][0] for x in repeat_types],
           [x.replace('_', ' ') for x in repeat_types],
           fontsize='small', handlelength=1, handleheight=1, handletextpad=0.5,
           frameon=False,
           bbox_to_anchor=(0,1.02,1,0.2), loc='lower left',
           mode='expand', borderaxespad=0, ncol=8)

# invert y-axis
ax.invert_yaxis()

sns.despine(left=True, offset=5, trim=True)

# save figure
fig = plt.gcf()

# without bbox_inches, the saved figure has truncated axes.
output_filename = 'stacked_barchart.meth_density.pdf'
fig.savefig(output_filename, bbox_inches='tight')

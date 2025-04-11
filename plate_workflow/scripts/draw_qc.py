import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

df = pd.read_csv(snakemake.input[0], index_col='wid')
order = [f'{row}{col:02}' for row in 'ABCDEFGHIJKLMNOP' for col in range(1,25)]
df = df.loc[order,:]

xs = list(range(1,25)) * 16
ys = [i for i in range(16,0,-1) for _ in range(24)]

for i,metric in enumerate(['total', 'uniq_frags_nuc', 'mapping_rate']):
    if metric == 'mapping_rate':
        colour = df[metric].values
    else:
        colour = np.log10(df[metric].values + 1)
    fig, axs = plt.subplots(figsize=(20,8), ncols=2, gridspec_kw={'width_ratios' : [13,7]})
    ax1, ax2 = axs
    sct = ax1.scatter(xs, ys, s=600, ec='k', c=colour)
    ax1.set_xticks(range(1,25))
    ax1.set_xticklabels([str(i) for i in range(1,25)], fontsize=20)
    ax1.set_yticks(range(1,17))
    ax1.set_yticklabels([i for i in 'ABCDEFGHIJKLMNOP'[::-1]], fontsize=20)
    ax1.set_xlim(0,25)
    ax1.set_ylim(0,17)
    ax1.tick_params(axis='x', labeltop=True, labelbottom=False, top=True, bottom=False)
    plt.colorbar(sct, ax=ax1, location='right', orientation='vertical', pad=0.01, fraction=0.1)

    ns, bins, ps = ax2.hist(colour, bins=50, density=True)

    fig.tight_layout()
    fig.savefig(snakemake.output[i], transparent=True, bbox_inches='tight')

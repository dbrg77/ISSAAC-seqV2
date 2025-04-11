import pandas as pd

summary = pd.read_csv(snakemake.input[0], index_col='barcode')
frip = pd.read_csv(snakemake.input[1], header=None,
                   names=['barcode', 'uniq_frags_in_peaks'], index_col='barcode')
chrm = pd.read_csv(snakemake.input[2], header=None,
                   names=['barcode', 'uniq_frags_mt'], index_col='barcode')

df = pd.concat([summary, frip, chrm], join='outer', axis=1)
df.fillna(0, inplace=True)

df['mapping_rate'] = round((df['total'] - df['unmapped'])/df['total'] * 100, 3)
df['uniq_frags_total'] = round(df['total'] - df['unmapped'] - df['duplicate'] - df['lowmapq'], 3)
df['uniq_frags_nuc'] = round(df['uniq_frags_total'] - df['uniq_frags_mt'], 3)
df['pct_mt'] = round(df['uniq_frags_mt']/df['uniq_frags_total'] * 100, 3)
df['frip'] = round(df['uniq_frags_in_peaks']/df['uniq_frags_total'], 3)

df.to_csv(snakemake.output[0], index=True)

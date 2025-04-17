import gzip
from scipy.io import mmwrite
from scipy.sparse import lil_matrix

bc_idx = {}
with open(snakemake.input[1]) as fh:
    for i,bc in enumerate(fh):
        bc_idx[bc.strip()] = i

num_bcs = len(bc_idx.keys())
num_pks = len(open(snakemake.input[0]).readlines())

mtx = lil_matrix((num_pks, num_bcs), dtype=int)

with gzip.open(snakemake.input[2], 'rt') as fh:
    for i,line in enumerate(fh):
        row_idx = i
        count_info = line.strip().split('\t')[1]
        items = count_info.split(',')
        col_idx = [bc_idx[bc_count.split(':')[0]] for bc_count in items]
        col_val = [int(bc_count.split(':')[1]) for bc_count in items]
        mtx[row_idx, col_idx] = col_val

mtx = mtx.tocsc()

mmwrite(snakemake.output[0], mtx, field='integer')

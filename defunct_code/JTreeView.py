from __future__ import print_function
import numpy as np
import pandas as pd
from scipy.cluster import hierarchy
from Utils import load_to_locals
import sys
from glob import glob

try:
    reload(sys)
except NameError:
    from importlib import reload

bad_cols = {'25u_emb4_sl10_FPKM', '25u_emb7_sl08_FPKM', '25u_emb7_sl09_FPKM'}

if 'gene_map' not in locals():
    gene_map = sorted(glob('prereqs/gene_map_table*.tsv'))[-1]
    fbgn_data = pd.read_table(gene_map, index_col=0,
                              na_values=['-', 'NaN', ''], keep_default_na=False,
                              skipfooter=3,
                                skiprows=5).dropna(how='all', axis=1)
    fbgn_lookup = dict(fbgn_data['primary_FBid'].dropna())
    fbgn_map = dict(fbgn_data['sequence_loc'].dropna())

def make_treeview_files(basename, data, clusters=None, do_cluster = False):
    if clusters is None and do_cluster:
        sys.stderr.write("Not yet implemented!\n")
        pass
    make_cdt_file(basename + '.cdt', data, clusters)
    make_gtr_file(basename + '.gtr', clusters)


def make_cdt_file(basename, data, clusters=None, sep_col = True):

    data = data.copy()
    if sep_col:
        prefixes = set(col[:col.find('_sl')] for col in data.columns)
        for prefix in prefixes:
            data[prefix+"_sep"] = pd.Series()
        data = data.sort_index(axis=1)

    data.insert(0, 'GID', 'NONE')
    data.insert(1, 'FBgn', data.index)
    data.insert(2, 'NAME', data.index)
    data.insert(3, 'CHROMOSOME', 'NONE')
    data.insert(4, 'ARM', 'L')
    data.insert(5, 'POSITION', 0)
    data.insert(6, 'GWEIGHT', 1.0)

    for i, row in enumerate(data.index):
        data.ix[row,'GID'] = 'GENE{}X'.format(i)
        data.ix[row, 'FBgn'] = fbgn_lookup.get(row, '???')
        if row in fbgn_map:
            pos = fbgn_map[row].split('..')[0]
            chrom, pos = pos.split(':')
            arm = 'R' if chrom.endswith('R') else 'L'
            if chrom[-1] in 'RL':
                chrom = chrom[:-1]
            data.ix[row, 'CHROMOSOME'] = chrom
            data.ix[row, 'ARM'] = arm
            data.ix[row, 'POSITION'] = int(pos)


    if clusters is not None:
        data = data.ix[hierarchy.leaves_list(clusters)]
    data.to_csv(basename, sep='\t', index=False, float_format='%.5f')

def make_gtr_file(basename, clusters):
    if isinstance(clusters, pd.DataFrame):
        pass
    else:
        clusters = pd.DataFrame(clusters,
                                columns=['left', 'right', 'score', 'children'])

    min = clusters.score.min()
    max = clusters.score.max()
    n = len(clusters) + 1
    def scale(score):
        return ((score - min)/(max - min+1e-10) * 2 - 1)*-1

    def format_name(idx):
        if idx < n:
            return 'GENE{:d}X'.format(int(idx))
        else:
            return 'NODE{:d}X'.format(int(idx - n + 1))

    datalist = []
    for i, row in clusters.iterrows():
        left = row['left']
        right = row['right']
        score = scale(row['score'])
        datalist.append(dict(NODEID='NODE{}X'.format(i+1),
                             LEFT=format_name(left),
                             RIGHT=format_name(right),
                             CORRELATION=score))

    out = pd.DataFrame(datalist, columns=['NODEID', 'LEFT', 'RIGHT', 'CORRELATION'])
    out.to_csv(basename, index=0, float_format='%.5f', sep='\t')
    return out

if __name__ == "__main__":
    import DistributionDifference
    reload(DistributionDifference)
    is_sparse = ''
    if '-sparse' in sys.argv:
        is_sparse='sparse_'
        try:
            step = int(sys.argv[sys.argv.index('-sparse') + 1])
        except:
            step = 10
        print('Sparse: {}'.format(step))
    else:
        step = 1

    expr_min = 15
    eps = 1.01
    read_table_args = dict(index_col=0,
                           keep_default_na=False,
                           na_values=['---', '', 'NA'])

    ase = pd.read_table('analysis_godot/ase_summary.tsv', **read_table_args).ix[::step]
    all_expr = pd.read_table('analysis_godot/summary_fb.tsv', **read_table_args).ix[::step]
    print("Read expression in")

    ase_rows, ase_cols = ase.shape
    ase = ase.ix[np.sum(np.isfinite(ase), axis=1) > .75 * ase_cols]
    ase = ase.ix[ase.index.intersection(all_expr.index)].dropna(axis='columns', how='all')
    all_expr = all_expr.ix[ase.index]

    all_expr_lognorm = np.log(all_expr+1).divide(np.log(all_expr.max( axis=1)+1),
                                                 axis=0)

    print("Precalculating distances")
    metric = DistributionDifference.earth_mover_multi
    dist_mat = DistributionDifference.mp_pandas_pdist(ase+eps, metric)

    Z = hierarchy.linkage(dist_mat, method='weighted')


    make_treeview_files("analysis/results/all_log_normed_"+is_sparse+metric.__name__,
                        all_expr_lognorm,
                        Z)

    make_treeview_files("analysis/results/all_"+is_sparse+metric.__name__,
                        all_expr,
                        Z)

    make_treeview_files("analysis/results/ase_"+is_sparse+metric.__name__,
            ase,
            Z)

    make_treeview_files("analysis/results/all_maxnorm_"+is_sparse+metric.__name__,
                        all_expr.divide(all_expr.max(axis=1) + 1, axis=0),
                        Z)


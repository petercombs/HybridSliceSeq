from __future__ import print_function
import pandas as pd
from scipy import stats
from collections import defaultdict

svase_dist_bed = 'analysis/results/transposon_dists/svase_dists.bed'
nosvase_dist_bed = 'analysis/results/transposon_dists/nosvase_dists.bed'
colnames = ['gchrom', 'gstart', 'gend', 'gname',
            'techrom', 'testart', 'teend',
            'dist']
if __name__ == "__main__":
    svase_dists = pd.read_table(svase_dist_bed,
                                header=None, names=colnames)
    nosvase_dists = pd.read_table(nosvase_dist_bed,
                                  header=None, names=colnames)

    closest_svase = defaultdict(lambda: 1e10)
    closest_nosvase = defaultdict(lambda: 1e10)

    for df, dist_dict in [(svase_dists, closest_svase),
                     (nosvase_dists, closest_nosvase)]:
        for i, row in df.iterrows():
            if row.dist == -1: continue
            old = dist_dict[row.gname]
            dist_dict[row.gname] = min(old, row.dist)

    sv_dist = pd.Series(index=closest_svase.keys(),data=pd.np.nan)
    nosv_dist = pd.Series(index=closest_svase.keys(),data=pd.np.nan)

    for gene in closest_nosvase:
        sv_gene = gene.split('<--')[1]
        sv_dist[sv_gene] = closest_svase[sv_gene]
        nosv_dist[sv_gene] = closest_nosvase[gene]

    print(stats.ttest_rel(sv_dist, nosv_dist, nan_policy='omit'))


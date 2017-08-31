import numpy as np
import pandas as pd
from matplotlib import pyplot as mpl
from collections import defaultdict
from Utils import sel_startswith, startswith, sel_contains
from scipy.stats import spearmanr
from progressbar import ProgressBar as pb


identity = lambda x: x
pearson = lambda x, y: np.corrcoef(x, y)[0,1]
spearman = lambda x, y: spearmanr(x, y)[0]
logp1 = lambda x: np.log10(x+1)

def get_corrs(data, adjust=identity, corr_func='pearson'):
    max_slice = defaultdict(int)
    for sl in data.columns:
        sl = sl.split('_sl')
        emb = sl[0]
        max_slice[emb] = max(max_slice[emb], int(sl[1][0:2]))
    xs = pd.Series(index=data.columns,
                   data=[int(a.split('_sl')[1][:2])/max_slice[a.split('_sl')[0]]
                         for a in data.columns if 'sl' in a])

    corrs_same = defaultdict(list)
    corrs_diff = defaultdict(list)
    all_corrs = [corrs_diff, corrs_same]
    for emb1_name in pb()(max_slice):
        emb1 = data.select(**sel_startswith(emb1_name)).applymap(adjust)
        genotype = emb1_name.split('_')[0]
        xs1 = xs.select(startswith(emb1_name))
        for emb2_name in max_slice:
            if emb1_name == emb2_name: continue
            emb2 = data.select(**sel_startswith(emb2_name)).applymap(adjust)
            xs2 = xs.select(startswith(emb2_name))
            closest = {
                column:
                min((abs(x2 - x1), c2)
                    for c2, x2 in xs2.items())[1]
                for column, x1 in xs1.items()
            }
            for col in emb1.columns:
                same = genotype == emb2_name.split('_')[0]
                all_corrs[same][genotype].append(emb1.ix[:, col].corr(
                    emb2.ix[:, closest[col]],
                    corr_func,
                ))
    return all_corrs





if __name__ == "__main__":
    ase = (pd
           .read_table('analysis_godot/ase_summary_by_read.tsv',
                       index_col=0,
                       keep_default_na=False, na_values=['---'],)
           .select(**sel_contains('X'))
           .dropna(how='all', axis=1)
           .dropna(how='all', axis=0)
          )

    expr = (pd
           .read_table('analysis_godot/summary.tsv',
                       index_col=0,
                       keep_default_na=False, na_values=['---'],)
           .dropna(how='all', axis=1)
           .dropna(how='all', axis=0)
          )

    all_ase_corrs = get_corrs(ase)
    all_expr_corrs = get_corrs(expr, adjust=np.log1p)


    #print(np.median(all_ase_corrs[0]), np.median(all_ase_corrs[1]))
    #mpl.hist(all_ase_corrs[0])

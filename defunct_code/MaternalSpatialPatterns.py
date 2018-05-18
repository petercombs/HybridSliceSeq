import DistributionDifference as dd
import pandas as pd
import Utils as ut
import PlotUtils as pu
from tqdm import tqdm
from warnings import filterwarnings

filterwarnings('ignore', category=FutureWarning)
filterwarnings('ignore', category=DeprecationWarning)

if __name__ == "__main__":
    expr = (pd
            .read_table('analysis_godot/summary.tsv', **ut.pd_kwargs)
            .drop('---', axis=1, errors='ignore')
           )
    ase = (pd
           .read_table('analysis_godot/wasp_summary_by_read.tsv',
                        **ut.pd_kwargs)
           .select(**ut.sel_startswith(('melXsim', 'simXmel')))
          )

    maternals = {line.strip() for line in open('analysis/results/maternal.txt')}

    mel = expr.select(**ut.sel_startswith('melXmel_'))
    sim = expr.select(**ut.sel_startswith('simXsim_'))
    hyb = expr.select(**ut.sel_startswith(('melXsim', 'simXmel')))

    emds = pd.Series(index=maternals, data=pd.np.nan)
    ones = pd.Series(index=expr.columns, data=1)
    for gene in tqdm(maternals):
        if gene not in expr.index: continue
        emds[gene] = dd.earth_mover_multi_rep(
            expr.ix[gene]+1, ones
        )
    emds = emds.sort_values(na_position='first')


    emds2 = emds[-50:]
    kwargs = pu.kwargs_expr_heatmap.copy()
    kwargs['progress_bar'] = True
    kwargs.pop('norm_rows_by')
    pu.svg_heatmap((expr.ix[emds2.index], ase.ix[emds2.index]),
                   'analysis/results/matspatpat.svg',
                   cmap=(pu.ISH, pu.cm.RdBu),
                   norm_rows_by=('maxall', 'center0pre'),
                   row_labels=[
                       ('{:.03f}'.format(emds[i]), i)
                       for i in emds2.index
                   ],
                   **kwargs)


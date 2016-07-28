"""Look for genes with different expression patterns in the hybrids

One thing to look for is genes where the two parents have fairly similar
expression patterns, but there is a noticeable shift in the hybrids.  This could
indicate some kind of cis + trans compensatory change.

"""

import numpy as np
import pandas as pd
import PlotUtils as pu
import DistributionDifference as dd
from matplotlib import cm
from progressbar import ProgressBar as pbar
from Utils import (pd_kwargs, sel_startswith)



EXPR_MIN = 5
NORMER = np.sum
pu_kwargs = {
    'box_height': 60,
    'col_sep': '_sl',
    'convert': True,
    'draw_box': True,
    'draw_name': False,
    'draw_row_labels': True,
    'make_hyperlinks': True,
    'max_width': 221,
    'progress_bar': False,
    'split_columns': True,
    'total_width': 200,
    'vspacer': 10}

def plot_expr_comparison(expr, gene, prefix=None, smoothed=0):
    mel = expr.select(**sel_startswith('mel_')).ix[gene]
    sim = expr.select(**sel_startswith('sim_')).ix[gene]
    hyb = expr.select(**sel_startswith(('melXsim', 'simXmel'))).ix[gene]

    if smoothed:
        mel = pd.rolling_mean(mel, smoothed, min_periods=1, center=True)
        sim = pd.rolling_mean(sim, smoothed, min_periods=1, center=True)
        hyb = pd.rolling_mean(hyb, smoothed, min_periods=1, center=True)
    pu.svg_heatmap(
        (None, mel, sim, hyb),
        'analysis_godot/results/spatial_diffs/{}.svg'.format(gene),
        cmap=(gene, cm.Reds,
               cm.Blues,
               pu.ISH),
        norm_rows_by=tuple([gene] + ['maxall']*7),
        **pu_kwargs
    )



if __name__ == "__main__":
    expr = pd.read_table('analysis_godot/summary.tsv', **pd_kwargs)
    ase = pd.read_table('analysis_godot/ase_summary_by_read.tsv', **pd_kwargs)

    mel = expr.select(**sel_startswith('mel_'))
    sim = expr.select(**sel_startswith('sim_'))
    hybrids = expr.select(**sel_startswith(('melXsim', 'simXmel')))
    melXsim = expr.select(**sel_startswith('melXsim'))
    simXmel = expr.select(**sel_startswith('simXmel'))

    expr_in_mel = (mel.max(axis=1) > EXPR_MIN)
    expr_in_sim = sim.max(axis=1) > EXPR_MIN
    expr_in_hybrids = (hybrids.max(axis=1) > EXPR_MIN)
    expr_in_all = (expr_in_mel & expr_in_sim & expr_in_hybrids)

    hyb_spatial_difference = pd.Series(data=np.nan,
                                       index=expr_in_all.index[expr_in_all]
                                      )
    parental_diffs = hyb_spatial_difference.copy()
    mel_hyb_diffs = hyb_spatial_difference.copy()
    sim_hyb_diffs = hyb_spatial_difference.copy()

    for gene in pbar()(hyb_spatial_difference.index):
        parental_diffs[gene] = dd.earth_mover_multi_rep(
            mel.ix[gene], sim.ix[gene],
            normer=lambda x: expr.ix[gene].max(),
        )
        mel_hyb_diffs[gene] = dd.earth_mover_multi_rep(
            mel.ix[gene], melXsim.ix[gene],
            normer=lambda x: expr.ix[gene].max(),
        )
        sim_hyb_diffs[gene] = dd.earth_mover_multi_rep(
            sim.ix[gene], simXmel.ix[gene],
            normer=lambda x: expr.ix[gene].max(),
        )
        hyb_spatial_difference[gene] = ((mel_hyb_diffs[gene] + sim_hyb_diffs[gene])
                                        / 2
                                        - parental_diffs[gene])



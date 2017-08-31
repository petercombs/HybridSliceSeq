"""Look for genes with different expression patterns in the hybrids

One thing to look for is genes where the two parents have fairly similar
expression patterns, but there is a noticeable shift in the hybrids.  This could
indicate some kind of cis + trans compensatory change.

"""

import numpy as np
import pandas as pd
import PlotUtils as pu
import HybridUtils as hu
import DistributionDifference as dd
from matplotlib import cm
from multiprocessing import Pool
from progressbar import ProgressBar as pbar
from Utils import (pd_kwargs, startswith, sel_startswith, get_synonyms, get_xs)
from CisTransASE import fit_all_splines



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


def get_diffs(expr, mel_spline, sim_spline, col_headers, offset=EXPR_MIN):
    mel = expr.select(startswith('melXmel_'))
    sim = expr.select(startswith('simXsim_'))
    melXsim = expr.select(startswith('melXsim_'))
    simXmel = expr.select(startswith('simXmel_'))
    hybrids = expr.select(startswith(('melXsim', 'simXmel')))
    parental_diffs = dd.earth_mover_multi_rep(
        mel+offset, sim+offset,
        #normer=lambda x: expr.max(),
    )
    mel_hyb_diffs = dd.earth_mover_multi_rep(
        mel+offset, melXsim+offset,
        #normer=lambda x: expr.max(),
    )
    sim_hyb_diffs = dd.earth_mover_multi_rep(
        sim+offset, simXmel+offset,
        #normer=lambda x: expr.max(),
    )

    hyb_hyb_diffs = dd.earth_mover_multi_rep(
        melXsim+offset, simXmel+offset,
        #normer=lambda x: expr.max(),
        #normer=pd.np.sum,
    )
    within_melXsim_diff = dd.earth_mover_within(
        melXsim+offset,
        #normer=expr.max(),
    )
    within_simXmel_diff = dd.earth_mover_within(
        simXmel+offset,
        #normer=expr.max(),
    )


    avgs = pd.Series((mel_spline(xs) + sim_spline(xs))/2,
                     index=col_headers,
                    )

    avg_hyb_diffs = dd.earth_mover_multi_rep(
        avgs.astype(float).clip(0, 1e6),
        hybrids,
        normer=lambda x: expr.max(),
    )
    avg_level = avgs.max()
    hyb_level = [hybrids.select(startswith(g)).max()
                 for g in ['melXsim_cyc14C_rep1', 'melXsim_cyc14C_rep2', 'melXsim_cyc14C_rep3',
                            'simXmel_cyc14C_rep1', 'simXmel_cyc14C_rep2']]

    return (
        hyb_hyb_diffs,
        parental_diffs, mel_hyb_diffs, sim_hyb_diffs, avgs, avg_hyb_diffs,
        avg_level, hyb_level,
        within_melXsim_diff, within_simXmel_diff,
    )

def plot_expr_comparison(expr, gene, prefix=None, smoothed=0):
    mel = expr.select(**sel_startswith('melXmel_')).ix[gene]
    sim = expr.select(**sel_startswith('simXsim_')).ix[gene]
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

def parse_args():
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--multi', default=True, action='store_true')
    parser.add_argument('--no-multi', dest='multi', action='store_false')
    return parser.parse_args()



if __name__ == "__main__":
    synonyms = get_synonyms()

    args = parse_args()
    expr = pd.read_table('analysis_godot/summary.tsv', **pd_kwargs).drop('---',
                                                                         axis='columns')
    ase = pd.read_table('analysis_godot/ase_summary_by_read.tsv',
                        **pd_kwargs).drop('---', axis='columns')


    mel = expr.select(**sel_startswith('melXmel_'))
    sim = expr.select(**sel_startswith('simXsim_'))
    hybrids = expr.select(**sel_startswith(('melXsim', 'simXmel')))
    melXsim = expr.select(**sel_startswith('melXsim'))
    simXmel = expr.select(**sel_startswith('simXmel'))

    expr_in_mel = (mel.max(axis=1) > EXPR_MIN)
    expr_in_sim = sim.max(axis=1) > EXPR_MIN
    expr_in_hybrids = (hybrids.max(axis=1) > EXPR_MIN)
    expr_in_all = (expr_in_mel & expr_in_sim & expr_in_hybrids)

    expr = expr.ix[expr_in_all]
    ase = ase.ix[expr.index]
    ase_classes = hu.get_classes(ase, pbar=pbar, style='cutoff')
    not_maternal = ase_classes.index[~((ase_classes.melXsim == 0) &
                                       (ase_classes.simXmel == 0))]

    with Pool() as p:
        mel_splines = fit_all_splines(mel.ix[not_maternal], p, progress=True)
        sim_splines = fit_all_splines(sim.ix[not_maternal], p, progress=True)

    hyb_spatial_difference = pd.Series(
        data=np.nan,
        index=not_maternal,
    )
    parental_diffs = hyb_spatial_difference.copy()
    hyb_hyb_diffs = hyb_spatial_difference.copy()
    mel_hyb_diffs = hyb_spatial_difference.copy()
    sim_hyb_diffs = hyb_spatial_difference.copy()
    avg_hyb_diffs = hyb_spatial_difference.copy()
    avg_levels = hyb_spatial_difference.copy()
    within_diffs_mXs = hyb_spatial_difference.copy()
    within_diffs_sXm = hyb_spatial_difference.copy()
    hyb_levels = pd.DataFrame(index=hyb_spatial_difference.index,
                           columns=['melXsim_cyc14C_rep1',
                                     'melXsim_cyc14C_rep2',
                                     'melXsim_cyc14C_rep3',
                                     'simXmel_cyc14C_rep1',
                                     'simXmel_cyc14C_rep2'])


    xs = np.linspace(0, 1, 20, endpoint=True)
    avgs = pd.DataFrame(index=hyb_spatial_difference.index,
                        columns=['avg_sl{}'.format(i+1) for i in range(20)])

    if args.multi:
        with Pool() as p:
            results = (p.apply_async(get_diffs, (expr.ix[gene], mel_splines[gene],
                                                 sim_splines[gene], avgs.columns))
                      for gene in hyb_spatial_difference.index)
            for gene in pbar()(hyb_spatial_difference.index):

                res = next(results).get()
                #res = get_diffs(expr.ix[gene], mel_splines[gene], sim_splines[gene],
                                #avgs.columns)
                (hyb_hyb_diffs[gene],
                 parental_diffs[gene],
                 mel_hyb_diffs[gene],
                 sim_hyb_diffs[gene],
                 avgs.ix[gene],
                 avg_hyb_diffs[gene],
                 avg_levels[gene], hyb_levels.ix[gene],
                 within_diffs_mXs[gene],
                 within_diffs_sXm[gene],
                ) = res
                hyb_spatial_difference[gene] = (avg_hyb_diffs[gene]
                                                - parental_diffs[gene])
    else:
        for gene in hyb_spatial_difference.index:
            print("About to try:", gene)
            res = get_diffs(expr.ix[gene], mel_splines[gene], sim_splines[gene],
                                avgs.columns)
            (hyb_hyb_diffs[gene],
             parental_diffs[gene],
             mel_hyb_diffs[gene],
             sim_hyb_diffs[gene],
             avgs.ix[gene],
             avg_hyb_diffs[gene],
             avg_levels[gene], hyb_levels.ix[gene],
             within_diffs_mXs[gene],
             within_diffs_sXm[gene],
            ) = res
            hyb_spatial_difference[gene] = (avg_hyb_diffs[gene]
                                            - parental_diffs[gene])
    l = locals()
    all_diffs = pd.DataFrame({i: l[i]
                              for i in ['hyb_hyb_diffs', 'parental_diffs',
                                        'mel_hyb_diffs', 'sim_hyb_diffs',
                                        'avg_hyb_diffs', 'within_diffs_mXs',
                                        'within_diffs_sXm',
                                       ]})



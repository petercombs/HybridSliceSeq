"""Determine genes that seem to be driven by cis and/or trans


If genes are mostly driven in cis and the trans environments are mostly the
same, then what we would expect to see is that the allele specific expression is
simply a function of the parental expression.

To test this, I'm going to fit a curve to the parental absolute expression data
separately (reported in FPKMs), then in the svASE samples, calculate the
difference between measured ASE and the predicted ASE, if it's just a function
of parental expression.

One particular thing that I may have to worry about here is the effect of
maternal expression.  I'm not sure as I start in on this what's the correct way
to deal with it, so I'm just going to ignore it until I start getting really
funny results.

"""
import numpy as np
import pandas as pd
from scipy import interpolate
from multiprocessing import Pool
from progressbar import ProgressBar as pbar
from Utils import pd_kwargs, sel_startswith, get_xs
from matplotlib import cm
import PlotUtils as pu
import DistributionDifference as dd


def fit_all_splines(expr, pool=None, progress=False):
    xs = get_xs(expr)
    is_good = (expr.isnull().sum() == 0)

    if pool is True:
        pool = Pool()

    asyncs = {gene: pool.apply_async(interpolate.UnivariateSpline,
                               (
                                   xs[is_good],
                                   pd.rolling_mean(
                                       expr.ix[gene, is_good],
                                       5,
                                       center=True,
                                       min_periods=1,
                                   ),
                                   None, [None, None],
                                   #3,
                                   #5,
                               )
                                    )
           for gene in expr.index}
    pb = pbar(maxval=len(asyncs) + 1)

    out = {}
    for i, gene in enumerate(asyncs):
        pb.update(i)
        res = asyncs[gene]
        out[gene] = res.get()
    pb.finish()
    return out


EXPR_MIN = 10
if __name__ == "__main__":
    print("Reading data")
    expr = (pd.read_table('analysis_godot/summary_fb.tsv', **pd_kwargs)
           .dropna(how='all', axis=0)
           )
    ase = (pd.read_table('analysis_godot/ase_summary_by_read.tsv', **pd_kwargs)
           .select(**sel_startswith(('melXsim', 'simXmel')))
           .dropna(how='all', axis=0)
          )
    paris = pd.read_table('prereqs/GSE68062_Gene_CLASS_after_FPKM_normalization.txt',
                  index_col=1)['mel.CLASS']
    pzyg = paris[paris != 'mat']
    mel = expr.select(**sel_startswith('mel_'))
    sim = expr.select(**sel_startswith('sim_'))
    hyb = expr.select(**sel_startswith(('melXsim', 'simXmel')))

    both_expr = (mel.max(axis=1) > EXPR_MIN) & (sim.max(axis=1) > EXPR_MIN)
    both_expr = (both_expr
                 .select(ase.index.__contains__)
                 .select(pzyg.index.__contains__)
                )
    both_expr = both_expr.index[both_expr]
    mel = mel.ix[both_expr]
    sim = sim.ix[both_expr]

    if 'mel_splines' not in locals() or locals().get('recalc', True):
        print("Fitting splines...")
        with Pool() as p:
            mel_splines = fit_all_splines(mel, p)
            sim_splines = fit_all_splines(sim, p)
        recalc = False

    ase_xs = get_xs(ase)
    ase_maternals = pd.Series(
        index=ase_xs.index,
        data=[1 if col.startswith('simXmel') else -1 for col in ase_xs.index]
    )
    ase_avgs = pd.DataFrame(
        data=dict(emd=np.nan, exprclass='?', avg_actual=np.nan,
                  avg_predicted=np.nan, bias=np.nan),
        index=mel.index,
    )

    prog = pbar(maxval=len(ase_avgs.index))
    #prog = lambda x: x
    for gene in prog(ase_avgs.index):
        good_ase = np.isfinite(ase.ix[gene]) & ~(ase.ix[gene] == ase_maternals)
        xg = ase_xs[good_ase]
        ase_avgs.ix[gene, 'actual'] = ase.ix[gene].mean()
        sim_pred = pd.Series(sim_splines[gene](ase_xs).clip(0.1, 1e10),
                             name='predicted_sim_'+gene,
                             index=ase_xs.index)
        mel_pred = pd.Series(mel_splines[gene](ase_xs).clip(0.1, 1e10),
                             name='predicted_mel_'+gene,
                             index=ase_xs.index)
        if sum(sim_pred[xg] + mel_pred[xg]) == 0:
            continue
        #print(gene)
        #print("Sim", min(sim_pred), max(sim_pred))
        #print("Mel", min(mel_pred), max(mel_pred))
        pred_ase = (
            ((sim_pred - mel_pred) /(sim_pred + mel_pred))
        )
        pred_ase.name='predicted_ase'
        '''
        print(pd.DataFrame([
            pd.Series(sim_pred, index=xg.index),
            pd.Series(mel_pred, index=xg.index),
            pred_ase,
            ase.ix[gene, good_ase],
        ],
            index=['sim_pred', 'mel_pred', 'ase_pred', 'ase_actual']
        ).T.to_csv(sep='\t', float_format='%.02f'))
        '''

        ase_avgs.ix[gene, 'predicted'] = np.average(
            pred_ase[good_ase],
            weights=(sim_pred + mel_pred)[good_ase],
        )
        ase_avgs.ix[gene, 'bias'] = (
            (ase.ix[gene, good_ase] - pred_ase[good_ase])
        ).mean()

        pu_kwargs = {
            'box_height': 60,
             'col_sep': '_sl',
             'convert': True,
             'draw_box': True,
             'draw_name': False,
             'draw_row_labels': True,
             'make_hyperlinks': True,
             'max_width': 1320,
             'progress_bar': False,
             'split_columns': True,
             'total_width': 200,
             'vspacer': 0}
        pu.svg_heatmap(
            (
                None, mel.ix[gene], None, None, None, None,
                None, mel_pred,
                None, sim.ix[gene], None, None, None, None,
                None, sim_pred,
                None, ase.ix[gene],
                None, pred_ase,
                None, hyb.ix[gene],
            ),
            'analysis_godot/results/ase_preds/{}.svg'.format(gene),
            norm_rows_by=('mel expression', 'max', '', '', '', '',
                          'predicted mel expression', 'max',
                          'sim expression', 'max', '', '', '', '',
                          'predicted sim expression', 'max',
                          'ASE', 'center0pre',
                          'predicted ASE', 'center0pre',
                          'hybrid expression', 'max',
                         ),
            cmap=(None, pu.ISH, None, None, None, None,
                  None, cm.Reds,
                  None, pu.ISH, None, None, None, None,
                  None, cm.Blues,
                  None, cm.RdBu,
                  None, cm.RdBu,
                  None, pu.ISH,
                 ),
            **pu_kwargs
        )
        '''
        if len(xg):
            ase_avgs.ix[gene, 'emd'] = dd.earth_mover_multi(
                ase.ix[gene, good_ase],
                pred_ase,
                #normer=lambda x: 1,
            )
            '''
        if gene in paris.index:
            ase_avgs.ix[gene, 'exprclass'] = paris[gene]


    ase_avgs.to_csv('analysis/results/ase_avgs.tsv', sep='\t')




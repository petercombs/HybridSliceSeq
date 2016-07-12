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
from scipy import interpolate, stats
from multiprocessing import Pool
from progressbar import ProgressBar as pbar
from Utils import pd_kwargs, sel_startswith, get_xs
from matplotlib import cm
import PlotUtils as pu
import DistributionDifference as dd

def wilson95_pref(ref, alt):
    """Lower bound of the 95% confidence interval

    Calculate the 95% confidence interval of the p, assuming a Bernoulli trial
    that gave the results REF and ALT.  Then, if that interval contains 50%,
    just use that, otherwise take the bound closer to 50%.  Finally, convert to
    a preference index [-1, 1], instead of a probability [0, 1] by multiplying
    by 2, and subtracting 1.

    See https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval.
    """

    z = 1.96
    n = ref + alt
    phat = alt/n

    plusminus = z * np.sqrt(1/n * phat * (1-phat) + 1/(4 * n**2) * z**2)

    p_plus = 1/(1+z**2/n) * (phat + z**2/(2*n) + plusminus)
    p_minus = 1/(1+z**2/n) * (phat + z**2/(2*n) - plusminus)

    if p_minus < 0.5 < p_plus:
        p = 0.5
    elif p_minus > 0.5:
        p = p_minus
    elif p_plus < 0.5:
        p = p_plus
    else:
        raise ValueError("I think I really screwed the pooch on this one")

    return 2 * p - 1

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
    'nan_replace' : 0.5,
    'vspacer': 10}


EXPR_MIN = 10
if __name__ == "__main__":
    print("Reading data")
    expr = (pd.read_table('analysis_godot/summary_fb.tsv', **pd_kwargs)
           .dropna(how='all', axis=0)
           )
    ase = (pd.read_table('analysis_godot/ase_summary_by_read.tsv', **pd_kwargs)
           .select(**sel_startswith(('melXsim', 'simXmel')))
           .dropna(how='all', axis=0)
           #.replace(pd.np.nan, 0)
          )

    chrom_of = {}
    for row in open('prereqs/gene_map_table_fb_2016_01.tsv'):
        if row.startswith('#') or not row.strip():
            continue
        data = row.split()
        chrom_of[data[1]] = data[-1].split(':')[0]

    males = ('melXsim_cyc14C_rep3', 'simXmel_cyc14C_rep2')
    on_x = [chrom_of[gene] == 'X' for gene in ase.index]
    is_male = [col.startswith(males) for col in ase.columns]
    ase.ix[on_x, is_male] = np.nan

    paris = pd.read_table('prereqs/GSE68062_Gene_CLASS_after_FPKM_normalization.txt',
                  index_col=1)['mel.CLASS']
    pzyg = paris[paris != 'mat']

    insitu_annots = pd.read_csv('prereqs/insitu_annot.csv', header=None,
                                names=['name', 'CG', 'FBgn', 'stage', 'annot'])
    ishzyg = (insitu_annots
              .query('stage == 2 and (annot in ("segmentally repeated", "gap", "pair rule"))')
              .FBgn.dropna().unique()
             )

    mel = expr.select(**sel_startswith('mel_'))
    sim = expr.select(**sel_startswith('sim_'))
    hyb = expr.select(**sel_startswith(('melXsim', 'simXmel')))

    ase_melXsim = ase.select(**sel_startswith('melXsim'))
    ase_simXmel = ase.select(**sel_startswith('simXmel'))

    frac_mXs_ase = pd.np.isfinite(ase_melXsim).sum(axis=1) / len(ase_melXsim.columns)
    frac_sXm_ase = pd.np.isfinite(ase_simXmel).sum(axis=1) / len(ase_simXmel.columns)

    ase_ks = pd.Series(index=ase.index, data=pd.np.nan)
    for gene in ase_ks.index:
        if (frac_mXs_ase[gene] >.2) and (frac_sXm_ase[gene] > .2):
            ase_ks[gene] = stats.ks_2samp(
                ase_melXsim.ix[gene].dropna(),
                ase_simXmel.ix[gene].dropna()).pvalue

    similar_ase = ase_ks.index[ase_ks > .05]

    ase_zyg = pd.Series(index=ase.index, data=np.nan)
    for gene in ase_zyg.index:
        if (frac_mXs_ase[gene] >.2) and (frac_sXm_ase[gene] > .2):
            tstat, pval = stats.ttest_ind(ase_melXsim.ix[gene].dropna(),
                                          ase_simXmel.ix[gene].dropna())
            tstat, pval2 = stats.ttest_ind(ase_melXsim.ix[gene].dropna(),
                                           -ase_simXmel.ix[gene].dropna())
            ase_zyg.ix[gene] = pval2/pval
    zyg_genes = ase_zyg.index[ase_zyg < 1e6]




    both_expr = (mel.max(axis=1) > EXPR_MIN) & (sim.max(axis=1) > EXPR_MIN)
    both_expr &= (frac_mXs_ase > .2) & (frac_sXm_ase > .2)
    both_expr = (both_expr
                 .select(ase.index.__contains__)
                 #.select(pzyg.index.__contains__)
                 .select(zyg_genes.__contains__)
                 #.select(similar_ase.__contains__)
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
        redraw = True

    ase_xs = get_xs(ase)
    ase_maternals = pd.Series(
        index=ase_xs.index,
        data=[1 if col.startswith('simXmel') else -1 for col in ase_xs.index]
    )
    ase_avgs = pd.DataFrame(
        data=dict(emd=np.nan, exprclass='?', actual=np.nan,
                  predicted=np.nan, bias=np.nan, n_good_slices=np.nan,
                  rmsdiff=np.nan),
        index=mel.index,
    )

    all_pred_ase_nan = pd.DataFrame(
        data=np.nan,
        index=ase.index,
        columns=ase.columns
    )

    all_mel_pred = all_pred_ase_nan.copy()
    all_sim_pred = all_pred_ase_nan.copy()

    prog = pbar(maxval=len(ase_avgs.index))
    #prog = lambda x: x
    for gene in prog(ase_avgs.index):
        if not locals().get('redraw', True):
            break
        good_ase = np.isfinite(ase.ix[gene]) & ~(ase.ix[gene] == ase_maternals)
        xg = ase_xs[good_ase]
        ase_avgs.ix[gene, 'n_good_slices'] = len(xg)
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
        #pred_ase = pd.Series(index=ase_xs.index, name='predicted_ase',
                             #data=pd.np.nan)
        #for ix in pred_ase.index:
            #pred_ase[ix] = wilson95_pref(10*mel_pred[ix], 10*sim_pred[ix])
        pred_ase = ((sim_pred - mel_pred) /(sim_pred + mel_pred))
        pred_ase.name='predicted_ase'
        pred_ase -= (pred_ase.mean() - ase.ix[gene].mean())
        pred_ase_nan = pred_ase.where(np.isfinite(ase.ix[gene]), np.nan)
        all_pred_ase_nan.ix[gene] = pred_ase_nan
        all_sim_pred.ix[gene] = sim_pred
        all_mel_pred.ix[gene] = mel_pred

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
        ase_avgs.ix[gene, 'rmsdiff'] = np.sqrt((ase.ix[gene, good_ase]- pred_ase[good_ase])
                                               .apply(np.square).mean())
        ase_avgs.ix[gene, 'bias'] = (
            (ase.ix[gene, good_ase] - pred_ase[good_ase])
        ).mean()

        ase_avgs.ix[gene, 'emd'] = dd.earth_mover_multi(
            (pd.rolling_mean(ase.ix[gene], 3, min_periods=1, center=True)
             / 2 + .5),
            pred_ase_nan/ 2 + .5,
            #normer=lambda x: 1,
            normer=pd.np.sum,
        )

        pu.svg_heatmap(
            (
                None, mel.ix[gene], None, None, None, None,
                None, mel_pred,
                None, sim.ix[gene], None, None, None, None,
                None, sim_pred,
                None, ase.ix[gene],
                None, pred_ase_nan,
                None, hyb.ix[gene],
            ),
            'analysis_godot/results/ase_preds/{}.svg'.format(gene),
            norm_rows_by=('mel expression', 'max', '', '', '', '',
                          'predicted mel expression', 'max',
                          'sim expression', 'max', '', '', '', '',
                          'predicted sim expression', 'max',
                          'ASE', 'center0pre',
                          'predicted ASE - EMD {:.03f}'.format(ase_avgs.ix[gene, 'emd']), 'center0pre',
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
        if gene in paris.index:
            ase_avgs.ix[gene, 'exprclass'] = paris[gene]


    redraw = False

    ase_avgs.to_csv('analysis/results/ase_avgs.tsv', sep='\t')
    trans = (ase_avgs.predicted - ase_avgs.actual) / 2**.5
    sorted_trans = abs(trans).sort_values()
    top_bottom = sorted_trans.index[:50] | sorted_trans.index[-50:]
    for gene in top_bottom:
        pu.svg_heatmap(
            (
                None, mel.ix[gene], None, None, None, None,
                None, all_mel_pred.ix[gene],
                None, sim.ix[gene], None, None, None, None,
                None, all_sim_pred.ix[gene],
                None, ase.ix[gene],
                None, all_pred_ase_nan.ix[gene],
                None, hyb.ix[gene],
            ),
            'analysis_godot/results/ase_preds/by_trans/{:+0.3f}-{}.svg'.format(trans.ix[gene], gene),
            norm_rows_by=('mel expression', 'max', '', '', '', '',
                          'predicted mel expression', 'max',
                          'sim expression', 'max', '', '', '', '',
                          'predicted sim expression', 'max',
                          'ASE', 'center0pre',
                          'predicted ASE - EMD {:.03f}'.format(ase_avgs.ix[gene, 'emd']), 'center0pre',
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

    rms_sorted = ase_avgs.rmsdiff.sort_values()
    n = len(rms_sorted) // 10

    for gene in (rms_sorted.index[:n] | rms_sorted.index[-n:]):
        pu.svg_heatmap(
            (
                None, mel.ix[gene], None, None, None, None,
                None, all_mel_pred.ix[gene],
                None, sim.ix[gene], None, None, None, None,
                None, all_sim_pred.ix[gene],
                None, ase.ix[gene],
                None, all_pred_ase_nan.ix[gene],
                None, hyb.ix[gene],
            ),
            'analysis_godot/results/ase_preds/by_rms/{:+0.3f}-{}.svg'.format(rms_sorted.ix[gene], gene),
            norm_rows_by=('mel expression', 'max', '', '', '', '',
                          'predicted mel expression', 'max',
                          'sim expression', 'max', '', '', '', '',
                          'predicted sim expression', 'max',
                          'ASE', 'center0pre',
                          'predicted ASE - EMD {:.03f}'.format(ase_avgs.ix[gene, 'emd']), 'center0pre',
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





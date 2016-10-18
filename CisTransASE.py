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
from Utils import (pd_kwargs, sel_startswith, get_xs, get_nearest_slice,
                   get_chroms, fbgns)
from matplotlib import cm
import PlotUtils as pu
import DistributionDifference as dd
import itertools as it
from collections import defaultdict

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

def rolling_weighted_average(values, weights, *args, **kwargs):
    weight_reshape = np.isfinite(values).multiply(weights[values.index
                                                          if len(values.shape) < 2
                                                          else values.columns])

    return (pd.rolling_sum(values.multiply(weight_reshape), *args, **kwargs)
            .divide(pd.rolling_sum(weight_reshape, *args, **kwargs)))

def fit_all_splines(expr, pool=None, progress=False):
    xs = get_xs(expr)
    is_good = (expr.isnull().sum() == 0)

    out = {}
    if progress:
        pb = pbar()
    else:
        pb = lambda x: x


    if pool is True:
        close = True
        pool = Pool()
    elif pool is None:
        for gene in pb(expr.index):
            expr_smooth = pd.rolling_mean(expr.ix[gene], 3, center=True,
                                          min_periods=1)
            is_good = ~expr_smooth.isnull()
            out[gene] = interpolate.UnivariateSpline(
                xs[is_good], expr_smooth[is_good]
            )
        return out
    else:
        close = False


    asyncs = {}
    for gene in expr.index:
        expr_smooth = pd.rolling_mean(expr.ix[gene], 3, center=True,
                                      min_periods=1)
        is_good = ~expr_smooth.isnull()
        asyncs[gene] = pool.apply_async(interpolate.UnivariateSpline, (
            xs[is_good],
            expr_smooth
        ))


    for gene in pb(asyncs):
        res = asyncs[gene]
        out[gene] = res.get()
    if close:
        pool.close()
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
    'nan_replace' : True,
    'vspacer': 10}

def calculate_spline_variance_explained(ase, splines, weights=None):
    xs = get_xs(ase)
    if weights is None:
        weights = np.ones_like(xs)
    if not callable(splines):
        var_obs = (((ase - splines) * weights)**2).sum()
        var_tot = (((ase - ase.mean()) * weights)**2).sum()
        return 1-(var_obs / var_tot)
    elif hasattr(splines, 'index'):
        r2 = pd.Series(index=splines.index, data=np.inf)
    elif hasattr(splines, 'keys'):
        r2 = pd.Series(index=list(splines.keys()), data=np.inf)
    else:
        return 1-(((ase - splines(xs))**2).sum() / ((ase - ase.mean())**2).sum())
    for ix in r2.index:
        r2.ix[ix] = 1-(
            ((ase.ix[ix] - splines[ix](xs))**2).sum()
            / ((ase.ix[ix] - ase.ix[ix].mean())**2).sum()
        )

    return r2

EXPR_MIN = 10
if __name__ == "__main__":
    print("Reading data")
    expr = (pd.read_table('analysis_godot/summary_wasp.tsv', **pd_kwargs)
            .dropna(how='all', axis=0)
            .dropna(how='all', axis=1)
           )
    if expr.index[0].startswith('FBgn'):
        expr.index = fbgns[expr.index]
    ase = (pd.read_table('analysis_godot/ase_summary_by_read_with_wasp.tsv', **pd_kwargs)
           .select(**sel_startswith(('melXsim', 'simXmel')))
           .dropna(how='all', axis=0)
           .rename_axis(lambda x: x.split('_ase')[0], axis=1)
           #.replace(pd.np.nan, 0)
          )
    ase = ase.ix[expr.index]

    read_counts = pd.read_table('analysis_godot/map_stats.tsv',
                                index_col='LongName')

    chrom_of = get_chroms()

    males = ('melXsim_cyc14C_rep3', 'simXmel_cyc14C_rep2')
    on_x = chrom_of[ase.index] == 'X'
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

    mel = expr.select(**sel_startswith(('melXmel_', 'mel_')))
    sim = expr.select(**sel_startswith(('simXsim_', 'sim_')))
    hyb = expr.select(**sel_startswith(('melXsim', 'simXmel')))

    ase_melXsim = ase.select(**sel_startswith('melXsim')).ix[mel.index]
    ase_simXmel = ase.select(**sel_startswith('simXmel')).ix[mel.index]

    frac_mXs_ase = ase_melXsim.count(axis=1) / len(ase_melXsim.columns)
    frac_sXm_ase = ase_simXmel.count(axis=1) / len(ase_simXmel.columns)

    ase_ks = pd.Series(index=ase.index.intersection(frac_mXs_ase.index), data=pd.np.nan)
    for gene in ase_ks.index:
        if (frac_mXs_ase[gene] >.2) and (frac_sXm_ase[gene] > .2):
            ase_ks[gene] = stats.ks_2samp(
                ase_melXsim.ix[gene].dropna(),
                ase_simXmel.ix[gene].dropna()).pvalue

    similar_ase = ase_ks.index[ase_ks > .05]

    ase_zyg = pd.Series(index=ase_ks.index, data=np.nan)
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
                 #.select(zyg_genes.__contains__)
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
                  r2=np.nan,
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

    all_pred_ase = []

    prog = pbar(maxval=len(ase_avgs.index))
    #prog = lambda x: x
    render_pool = Pool()
    renders = []
    for gene in prog(ase_avgs.index):
        if not locals().get('redraw', True):
            break
        good_ase = np.isfinite(ase.ix[gene]) & ~(ase.ix[gene] == ase_maternals)
        xg = ase_xs[good_ase]
        ase_avgs.ix[gene, 'n_good_slices'] = len(xg)
        ase_avgs.ix[gene, 'actual'] = ase.ix[gene].mean()
        sim_pred = pd.Series(sim_splines[gene](ase_xs).clip(1e-3, 1e10),
                             name='predicted_sim_'+gene,
                             index=ase_xs.index)
        mel_pred = pd.Series(mel_splines[gene](ase_xs).clip(1e-3, 1e10),
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
        min_ase = ase.ix[gene].min()
        max_ase = ase.ix[gene].max()
        pred_ase = ((pred_ase - pred_ase.mean() + ase.ix[gene].mean())
                    * (max_ase - min_ase)/(pred_ase.max() - pred_ase.min()))
        pred_ase.name='predicted_ase'
        #pred_ase -= (pred_ase.mean() - ase.ix[gene].mean())
        pred_ase_nan = pred_ase.where(np.isfinite(ase.ix[gene]), np.nan)
        all_pred_ase_nan.ix[gene] = pred_ase_nan
        all_sim_pred.ix[gene] = sim_pred
        all_mel_pred.ix[gene] = mel_pred
        all_pred_ase.extend(
            (a, b) for  a, b in zip(
                ase.ix[gene, good_ase],
                pred_ase_nan[good_ase])
            if (np.isfinite(a) and np.isfinite(b)))

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

        ase_avgs.ix[gene, 'r2'] = calculate_spline_variance_explained(
            ase.ix[gene],
            pred_ase_nan
        )
        '''
        ase_avgs.ix[gene, 'emd'] = dd.earth_mover_multi(
            (pd.rolling_mean(ase.ix[gene], 3, min_periods=1, center=True)
             / 2 + .5),
            pred_ase_nan/ 2 + .5,
            #normer=lambda x: 1,
            normer=pd.np.sum,
        )
        '''

        renders.append(render_pool.apply_async(
            pu.svg_heatmap,
            (
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
            ),
            dict(
                norm_rows_by=('mel expression', 'max', '', '', '', '',
                              'predicted mel expression', 'max',
                              'sim expression', 'max', '', '', '', '',
                              'predicted sim expression', 'max',
                              'ASE', 'center0pre',
                              'predicted ASE - r2 {:.03f}'.format(ase_avgs.ix[gene, 'r2']), 'center0pre',
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
            ),
        ))
        if gene in paris.index:
            ase_avgs.ix[gene, 'exprclass'] = paris[gene]


    print(ase_avgs.ix[:, 'r2'].sort_values())
    redraw = False

    recalc_cis_baseline = locals().get('recalc_cis_baseline', True)

    embryos = {col.split('_sl')[0] for col in ase.columns}
    combos = list(it.permutations(embryos, 2))
    smoothed_ase_vals = defaultdict(list)
    actual_ase_vals = defaultdict(list)
    by_nearest = []
    for emb1, emb2 in pbar()(combos):
        if not recalc_cis_baseline:
            print("Skipping cis comparisons")
            break
        emb1 = ase.select(**sel_startswith(emb1))
        emb2 = ase.select(**sel_startswith(emb2))
        emb1_xs = get_xs(emb1)
        emb2_xs = get_xs(emb2)
        for slice1 in emb1.columns:
            slice2 = get_nearest_slice(emb1_xs[slice1], emb2)
            by_nearest.extend(
                (a, b) for gene, a, b in zip(
                    emb1.index,
                    emb1.ix[:, slice1],
                    emb2.ix[:, slice2])
                if (np.isfinite(a) and np.isfinite(b) and gene in both_expr))
        for gene in emb1.index:
            if (
                gene not in both_expr or
                sum(np.isfinite(emb1.ix[gene])) < 5 or
                sum(np.isfinite(emb2.ix[gene])) < 5 or
                False
               ):
                continue
            yvals = pd.rolling_mean( emb1.ix[gene], 3, center=True, min_periods=1)
            gx = np.isfinite(yvals)
            spline = interpolate.UnivariateSpline(
                emb1_xs[gx],
                yvals[gx],
                bbox=[0,1],
            )

            gx2 = (
                np.isfinite(emb2.ix[gene])
                & (emb2_xs > emb1_xs[gx][0])
                & (emb2_xs < emb1_xs[gx][-1])
                  )
            smoothed_vals = spline(emb2_xs[gx2])
            svnm = np.nanmean(smoothed_vals)
            if np.isnan(svnm):
                print(emb1.columns[0].split('_sl')[0],
                      emb2.columns[0].split('_sl')[0],
                      gene)
                #assert False
            smoothed_ase_vals[gene].extend(smoothed_vals)

            actual_ase_vals[gene].extend(emb2.ix[gene,gx2])



    recalc_cis_baseline = False

    ase_avgs.to_csv('analysis/results/ase_avgs.tsv', sep='\t')
    trans = (ase_avgs.predicted - ase_avgs.actual) / 2**.5
    redraw_outliers = locals().get('redraw_outliers', True)
    sorted_trans = abs(trans).sort_values()
    top_bottom = sorted_trans.index[:50] | sorted_trans.index[-50:]
    for gene in pbar()(top_bottom):
        if not redraw_outliers:
            break
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

    for gene in pbar()(rms_sorted.index[:n] | rms_sorted.index[-n:]):
        if not redraw_outliers:
            break
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
        redraw_outliers = False


    co = 1
    trans_by_gene = defaultdict(list)
    for gene in smoothed_ase_vals:
        for smo, act in zip(smoothed_ase_vals[gene], actual_ase_vals[gene]):
            trans_by_gene[gene].append((smo - act)/2**.5)
    frac_very_trans = {gene:
                       sum(abs(np.array(trans_by_gene[gene])) > .5)
                       / len(trans_by_gene[gene])
                       for gene in trans_by_gene}


    for item in pbar()(renders):
        item.get()
    render_pool.close()





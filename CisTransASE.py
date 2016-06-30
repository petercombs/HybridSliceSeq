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



def fit_all_splines(expr, pool=None, progress=False):
    xs = get_xs(expr)
    is_good = (expr.isnull().sum() == 0)

    if pool is True:
        pool = Pool()

    asyncs = {gene: pool.apply_async(interpolate.interp1d,
                               (xs[is_good],
                                expr.ix[gene, is_good],
                                "cubic"))
           for gene in expr.index}
    pb = pbar(maxval=len(asyncs) + 1)

    out = {}
    for i, gene in enumerate(asyncs):
        pb.update(i)
        res = asyncs[gene]
        out[gene] = res.get()
    pb.finish()
    return out


EXPR_MIN = 5
if __name__ == "__main__":
    print("Reading data")
    expr = (pd.read_table('analysis_godot/summary_fb.tsv', **pd_kwargs)
           .dropna(how='all', axis=0)
           .dropna(how='all', axis=1))
    ase = (pd.read_table('analysis_godot/ase_summary_by_read.tsv', **pd_kwargs)
           .dropna(how='all', axis=0)
           .dropna(how='all', axis=1))
    paris = pd.read_table('prereqs/GSE68062_Gene_CLASS_after_FPKM_normalization.txt',
                  index_col=1)['mel.CLASS']
    mel = expr.select(**sel_startswith('mel_'))
    sim = expr.select(**sel_startswith('sim_'))

    both_expr = (mel.max(axis=1) > EXPR_MIN) & (sim.max(axis=1) > EXPR_MIN)
    both_expr = both_expr.select(ase.index.__contains__)
    both_expr = both_expr.index[both_expr]
    mel = mel.ix[both_expr]
    sim = sim.ix[both_expr]
    print("Fitting splines...")

    with Pool() as p:
        mel_splines = fit_all_splines(mel, p)
        sim_splines = fit_all_splines(sim, p)

    ase_xs = get_xs(ase)
    ase_avgs = pd.DataFrame(
        data=dict(exprclass='?', actual=np.nan, predicted=np.nan),
        index=mel.index,
    )

    for gene in pbar()(ase_avgs.index):
        xg = ase_xs[np.isfinite(ase.ix[gene])]
        ase_avgs.ix[gene, 'actual'] = ase.ix[gene].mean()
        sim_pred = sim_splines[gene](xg)
        mel_pred = mel_splines[gene](xg)
        ase_avgs.ix[gene, 'predicted'] = np.average(
            ((sim_pred - mel_pred) /(sim_pred + mel_pred)),
            weights=(sim_pred + mel_pred),
        )
        if gene in paris.index:
            ase_avgs.ix[gene, 'exprclass'] = paris[gene]


    ase_avgs.to_csv('analysis/results/ase_avgs.tsv', sep='\t')




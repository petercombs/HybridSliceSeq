from numpy import (exp, sqrt, array, isfinite, mean, shape, nan, inf, median)
from scipy.optimize import curve_fit
from scipy.stats import linregress
from multiprocessing import Pool
from collections import defaultdict
import pandas as pd
import itertools as it
import numpy as np
import inspect
from progressbar import ProgressBar

def logistic(x, A, k, x0, y0):
    return A/(1+exp(k*(x-x0))) - y0

def peak(x, A, w, x0, y0):
    return A * exp(-(x-x0)**2/w**2) - y0

def fit_func(func, index, data, xs, p0=None, randomize=False):
    ys = data.ix[index]
    if randomize:
        xs = np.random.permutation(xs)
    keep = array(isfinite(ys))
    if sum(keep) < len(xs)/2:
        return array([nan, nan, nan, nan])
    if p0 is None:
        pos_amp = ys.max() - ys.median()
        neg_amp = ys.min() - ys.median()
        if pos_amp > -neg_amp:
            amp = pos_amp
        else:
            amp = neg_amp
        p0 = [amp, 20, xs.mean(), ys.median()]
    if func == peak:
        w_min = 0.15
        w_max = 0.6
        p0[1] = 0.4
        if (ys.max() - ys.median()) < (ys.median() - ys.min()):
            p0[0] = ys.median() - ys.min()
            p0[2] = xs[ys.argmin()]
        else:
            p0[2] = xs[ys.argmax()]
    else:
        w_min = 1
        w_max = inf
    try:
        return curve_fit(func,
                         xs[keep],
                         ys[keep],
                         p0,
                         bounds=[[-2, w_min, -0.1, -1],
                                 [2, w_max, 1.1, 1]],
                        )[0]
    except (ValueError, RuntimeError):
        return array([nan, nan, nan, nan])

def fit_all_ase(ase, func, xs, colnames=None, pool=None):
    if colnames is None:
        colnames = inspect.getargspec(func).args
        colnames.pop('x')
    if pool is not None:
        res = list(pool.starmap(
            fit_func,
            zip(
                it.repeat(func),
                ase.index,
                it.repeat(ase),
                it.repeat(xs)
                )
            ))
    else:
        res = list(it.starmap(
            fit_func,
            zip(
                it.repeat(func),
                ase.index,
                it.repeat(ase),
                it.repeat(xs)
                )
            ))

    return pd.DataFrame(
            index=ase.index,
            columns=colnames,
            data=res
            )


def calculate_variance_explained(ase, xs, func, params):
    r2 = pd.Series(index=params.index, data=np.inf)
    for ix in params.index:
        r2.ix[ix] = 1-(
            ((ase.ix[ix] - func(xs, *params.ix[ix]))**2).sum()
            / ((ase.ix[ix] - ase.ix[ix].mean())**2).sum()
        )

    return r2


if __name__ == "__main__":
    ase = (pd
           .read_table('analysis_godot/ase_summary_by_read.tsv',
                       index_col=0,
                       keep_default_na=False, na_values=['---'],)
           .dropna(how='all', axis=1)
           .dropna(how='all', axis=0)
          )
    ase_perm = pd.DataFrame(
            data=np.random.permutation(ase.T).T,
            index=ase.index,
            columns=ase.columns,
            )
    max_slice = defaultdict(int)
    for sl in ase.columns:
        sl = sl.split('_')
        emb = '_'.join(sl[:2])
        max_slice[emb] = max(max_slice[emb], int(sl[2][2:]))

    xs = pd.Series(index=ase.columns,
                   data=[int(a.split('_')[2][2:])/max_slice['_'.join(a.split('_')[:2])] for a in ase.columns if 'sl' in a])
    colnames = ['Amp', 'width', 'center', 'y_offset']
    with Pool() as p:
        res_logist = fit_all_ase(ase, logistic, xs, colnames, p).dropna()
        res_logist_perm = fit_all_ase( ase_perm, logistic, xs, colnames, p).dropna()
        res_peak = fit_all_ase(ase, peak, xs, colnames, p).dropna()
        res_peak_perm = fit_all_ase( ase_perm, peak, xs, colnames, p).dropna()

        res_lin = pd.DataFrame(
                index=ase.index,
                columns=['slope', 'intercept', 'pval', 'r'],
                data=np.nan
                )
        res_lin_perm = res_lin.copy()

        for gene in ProgressBar()(ase.index):
            cols = isfinite(ase.ix[gene])
            linreg = linregress(xs[cols], ase.ix[gene, cols])
            if linreg.pvalue < .05:
                res_lin.ix[gene] = [linreg.slope, linreg.intercept, linreg.pvalue, linreg.rvalue]

            cols = isfinite(ase_perm.ix[gene])
            linreg = linregress(xs[cols], ase_perm.ix[gene, cols])
            if linreg.pvalue < .05:
                res_lin_perm.ix[gene] = [linreg.slope, linreg.intercept, linreg.pvalue, linreg.rvalue]



    r2_logist = calculate_variance_explained(ase, xs, logistic, res_logist)
    r2_logist_perm = calculate_variance_explained(ase_perm, xs, logistic, res_logist_perm)
    r2_peak = calculate_variance_explained(ase, xs, peak, res_peak)
    r2_peak_perm = calculate_variance_explained(ase_perm, xs, peak, res_peak_perm)


    print("Found {: 4} {:<10} of which {} are likely false".format(sum(r2_logist > 0.5), 'logistic', sum(r2_logist_perm > .5)))
    print("Found {: 4} {:<10} of which {} are likely false".format(sum(r2_peak > 0.5), 'peak', sum(r2_peak_perm > .5)))

    good_amps_logist = res_logist.Amp[r2_logist>0.5].sort_values(inplace=False)
    good_amps_peak = res_peak.Amp[r2_peak>0.5].sort_values(inplace=False)
    #good_amps_logist.to_csv('analysis/results/logist.tsv', sep='\t')
    #good_amps_peak.to_csv('analysis/results/peak.tsv', sep='\t')
    res_logist.ix[r2_logist > 0.5].to_csv('analysis/results/logist.tsv', sep='\t')
    res_peak.ix[r2_peak > 0.5].to_csv('analysis/results/peak.tsv', sep='\t')






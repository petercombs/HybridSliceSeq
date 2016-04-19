from numpy import (exp, sqrt, array, isfinite, mean, shape, nan, inf, median)
from scipy.optimize import curve_fit
from multiprocessing import Pool
from collections import defaultdict
import pandas as pd
import itertools as it
import numpy as np

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
        p0 = [ys.max() - ys.median(), 0.5, xs.mean(), ys.median()]
    if func == peak:
        w_min = 0
        w_max = 0.6
        if (ys.max() - ys.median()) < (ys.median() - ys.min()):
            p0[0] = ys.median() - ys.min()
            p0[2] = xs[ys.argmin()]
        else:
            p0[2] = xs[ys.argmax()]
    else:
        w_min = -inf
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

if __name__ == "__main__":
    ase = (pd
           .read_table('analysis_godot/ase_summary.tsv',
                       index_col=0,
                       keep_default_na=False, na_values=['---'],)
           .dropna(how='all', axis=1)
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
    with Pool() as p:
        res_logist_old = list(p.starmap(fit_func,
                               zip(it.repeat(logistic),
                                   ase.index,
                                   it.repeat(ase),
                                   it.repeat(xs),
                                  )
                              ))
        res_logist = pd.DataFrame(
            index=ase.index,
            columns=['Amp', 'width', 'center', 'y_offset'],
            data=res_logist_old,
        ).dropna()
        res_logist_perm = list(p.starmap(fit_func,
                               zip(it.repeat(logistic),
                                   ase_perm.index,
                                   it.repeat(ase_perm),
                                   it.repeat(xs),
                                  )
                              ))
        res_logist_perm = pd.DataFrame(
            index=ase_perm.index,
            columns=['Amp', 'width', 'center', 'y_offset'],
            data=res_logist_perm,
        ).dropna()
        res_peak = list(p.starmap(fit_func,
                               zip(it.repeat(peak),
                                   ase.index,
                                   it.repeat(ase),
                                   it.repeat(xs),
                                  )
                              ))
        res_peak =pd.DataFrame(
            index=ase.index,
            columns=['Amp', 'width', 'center', 'y_offset'],
            data=list(res_peak),
        ).dropna()
        res_peak_perm = list(p.starmap(fit_func,
                               zip(it.repeat(peak),
                                   ase_perm.index,
                                   it.repeat(ase_perm),
                                   it.repeat(xs),
                                  )
                              ))
        res_peak_perm =pd.DataFrame(
            index=ase_perm.index,
            columns=['Amp', 'width', 'center', 'y_offset'],
            data=list(res_peak_perm),
        ).dropna()

    r2_logist = pd.Series(index=res_logist.index, data=np.inf)
    r2_logist_perm = pd.Series(index=res_logist_perm.index, data=np.inf)
    r2_peak = pd.Series(index=res_peak.index, data=np.inf)
    r2_peak_perm = pd.Series(index=res_peak_perm.index, data=np.inf)

    for ix in r2_logist.index:
        r2_logist.ix[ix] = 1-(
            ((ase.ix[ix] - logistic(xs, *res_logist.ix[ix]))**2).sum()
            / ((ase.ix[ix] - ase.ix[ix].mean())**2).sum()
        )

    for ix in r2_peak.index:
        r2_peak.ix[ix] = 1-(
            ((ase.ix[ix] -     peak(xs,   *res_peak.ix[ix]))**2).sum()
            / ((ase.ix[ix] - ase.ix[ix].mean())**2).sum()
        )

    for ix in r2_logist_perm.index:
        r2_logist_perm.ix[ix] = 1-(
            ((ase_perm.ix[ix] - logistic(xs, *res_logist_perm.ix[ix]))**2).sum()
            / ((ase_perm.ix[ix] - ase_perm.ix[ix].mean())**2).sum()
        )

    for ix in r2_peak_perm.index:
        r2_peak_perm.ix[ix] = 1-(
            ((ase_perm.ix[ix] -     peak(xs,   *res_peak_perm.ix[ix]))**2).sum()
            / ((ase_perm.ix[ix] - ase_perm.ix[ix].mean())**2).sum()
        )

    print("Found {: 4} {:<10} of which {} are likely false".format(sum(r2_logist > 0.5), 'logistic', sum(r2_logist_perm > .5)))
    print("Found {: 4} {:<10} of which {} are likely false".format(sum(r2_peak > 0.5), 'peak', sum(r2_peak_perm > .5)))

    good_amps_logist = res_logist.Amp[r2_logist>0.5].sort_values(inplace=False)
    good_amps_peak = res_peak.Amp[r2_peak>0.5].sort_values(inplace=False)
    #good_amps_logist.to_csv('analysis/results/logist.tsv', sep='\t')
    #good_amps_peak.to_csv('analysis/results/peak.tsv', sep='\t')
    res_logist.ix[r2_logist > 0.5].to_csv('analysis/results/logist.tsv', sep='\t')
    res_peak.ix[r2_peak > 0.5].to_csv('analysis/results/peak.tsv', sep='\t')






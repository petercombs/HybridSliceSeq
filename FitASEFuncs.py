from numpy import (exp, sqrt, array, isfinite, mean, shape, nan, inf)
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

def fit_func(func, index, data, xs):
    ys = data.ix[index]
    keep = array(isfinite(ys))
    if sum(keep) < len(xs)/2:
        return array([nan, nan, nan, nan])
    p0 = [0.5, 1, mean(xs), mean(ys)]
    assert shape(xs[keep]) == shape(ys[keep])
    try:
        return curve_fit(func,
                         xs[keep],
                         ys[keep],
                         p0,
                         bounds=[[-1, -inf, 0.1, -1],
                                 [1, inf, 0.9, 1]],
                        )[0]
    except:
        return array([nan, nan, nan, nan])

if __name__ == "__main__":
    ase = (pd
           .read_table('ase_summary.tsv',
                       index_col=0,
                       keep_default_na=False, na_values=['---'],)
           .dropna(how='all', axis=1)
          )
    max_slice = defaultdict(int)
    for sl in ase.columns:
        sl = sl.split('_')
        emb = '_'.join(sl[:2])
        max_slice[emb] = max(max_slice[emb], int(sl[2][2:]))

    xs = array([int(a.split('_')[2][2:])/max_slice['_'.join(a.split('_')[:2])] for a in ase.columns if 'sl' in a])
    with Pool() as p:
        res_logist = list(p.starmap(fit_func,
                               zip(it.repeat(logistic),
                                   ase.index,
                                   it.repeat(ase),
                                   it.repeat(xs),
                                  )
                              ))
        res_logist = pd.DataFrame(
            index=ase.index,
            columns=['Amp', 'width', 'center', 'y_offset'],
            data=res_logist,
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

    r2_logist = pd.Series(index=res_logist.index, data=np.inf)
    r2_peak = pd.Series(index=res_peak.index, data=np.inf)

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

    good_amps_logist = res_logist.Amp[r2_logist>0.5].sort_values(inplace=False)






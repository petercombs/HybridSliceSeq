import numpy as np
from scipy import interpolate

import collections
import functools
import emd
import math
from Utils import contains
import pandas as pd
from itertools import combinations

class memoized(object):
   '''Decorator. Caches a function's return value each time it is called.
   If called later with the same arguments, the cached value is returned
   (not reevaluated).
   '''
   def __init__(self, func):
      self.func = func
      self.cache = {}
   def __call__(self, *args):
      if not isinstance(args, collections.Hashable):
         # uncacheable. a list, for instance.
         # better to not cache than blow up.
         return self.func(*args)
      if tuple(*args) in self.cache:
         return self.cache[tuple(*args)]
      else:
         value = self.func(*args)
         self.cache[tuple(*args)] = value
         return value
   def __repr__(self):
      '''Return the function's docstring.'''
      return self.func.__doc__
   def __get__(self, obj, objtype):
      '''Support instance methods.'''
      return functools.partial(self.__call__, obj)

@memoized
def convert_to_distribution(points):
    """Converts the given data points to a smoothed distribution from 0-100%

    """

    x = np.linspace(0,1, len(points), endpoint=True)
    f = interpolate.interp1d(x, points, kind='cubic')
    retval =  np.cumsum(f(np.linspace(0, 1, 30, endpoint=True)).clip(0,1e5))
    return retval / (sum(retval)+1e-10)

def get_nonuniform_mp(args, eps=1.01):
    gene, row = args
    temp = pd.Series(data=1, index=row.index)
    if sum(pd.np.isfinite(row))==0:
        return pd.np.nan
    return earth_mover_multi(row + eps, temp)

def get_nonuniform_inv_mp(args):
    gene, row = args
    temp = pd.Series(data=1, index=row.index)
    if sum(pd.np.isfinite(row))==0:
        return pd.np.nan
    return earth_mover_multi(eps - row, temp)

def diff_stat(points1, points2):

    dist1 = convert_to_distribution(points1)
    dist2 = convert_to_distribution(points2)


    normfac = np.log(max(max(points1), max(points2)) + 1)

    return np.max(np.abs(dist1 - dist2)) * normfac

divmat = np.zeros([0,0])


def tang_stat(points1, points2):
    assert len(points1) == len(points2)
    points1 = np.array(points1 / np.mean(points1))
    points2 = np.array(points2 / np.mean(points2))

    va = np.reshape(np.repeat(points1, len(points2)), (len(points2), -1),
                    order='C')
    vb = np.reshape(np.repeat(points2, len(points1)), (-1, len(points2)),
                    order='F')

    global divmat
    if np.shape(divmat) != (len(points1), len(points2)):
        x, y = np.mgrid[0:len(points1), 0:len(points2)]
        divmat = 1/(np.abs(x - y) + 1)
    return np.sqrt(np.sum(np.triu((va - vb)**2 * divmat)))

#    stat = 0
#    for i in range(len(points1)):
#        for j in range(len(points2)):
#            stat += (points1[i] - points2[j])**2 / (np.abs(i - j)+1)
#            if i == j:
#                break
#
#    return np.sqrt(stat)

def earth_mover(points1, points2, normer=np.sum):
    xs1 = np.linspace(0,1,len(points1),
                      endpoint=True)[np.array(np.isfinite(points1))]
    xs2 = np.linspace(0,1,len(points2),
                      endpoint=True)[np.array(np.isfinite(points2))]
    points1 = points1[np.isfinite(points1)]
    points2 = points2[np.isfinite(points2)]
    return emd.emd(xs1, xs2,
                   points1/normer(points1),
                   points2/normer(points2))

def lcm(a,b):
    return abs(a * b) / math.gcd(a,b) if a and b else 0

def earth_mover_interp(points1, points2, normer=np.sum):
    xs1 = np.linspace(0,1,len(points1),
                      endpoint=True)[np.array(np.isfinite(points1))]
    xs2 = np.linspace(0,1,len(points2),
                      endpoint=True)[np.array(np.isfinite(points2))]
    xs = np.linspace(0, 1, min(100, lcm(len(points1),len(points2))),
                     endpoint=True)
    if np.sum(np.isfinite(points1)) == 0:
        return 0
    points1 = np.interp(xs, xs1, points1[np.isfinite(points1)])
    points2 = np.interp(xs, xs2, points2[np.isfinite(points2)])
    return emd.emd(xs, xs,
                   points1/(normer(points1) if callable(normer) else normer),
                   points2/(normer(points2) if callable(normer) else normer))


startswith = lambda x: lambda y: y.startswith(x)

def earth_mover_multi_rep(points1, points2, normer=np.sum):
    dist = 0.0
    reps1 = {col.split('sl')[0] for col in
             points1.index
            }
    reps2 = {col.split('sl')[0] for col in
             points2.index
            }
    for rep1 in reps1:
        for rep2 in reps2:
            dist += (earth_mover_interp(points1.select(contains(rep1)),
                                        points2.select(contains(rep2)),
                                        normer=normer)**2
                     / (len(reps1)*len(reps2)))
    return dist**.5

def earth_mover_within(points, sep='_sl', normer=np.sum):
   dist = 0.0
   n = 0
   reps = {
       col.split(sep)[0] for col in points.index
   }
   for points1, points2 in combinations(reps, 2):
       n += 1
       dist += earth_mover_interp(points.select(contains(points1)),
                                   points.select(contains(points2)),
                                   normer=normer)
   return dist / n


def earth_mover_multi(points1, points2, normer=np.sum):
    dist = 0.0
    embs = {col.split('sl')[0] for col in points1.index}
    #sums = [[],[]]
    for emb in embs:
        if ((sum(np.isfinite(points1.select(startswith(emb)))) == 0)
            or (sum(np.isfinite(points2.select(startswith(emb)))) == 0)):
            continue
        dist += earth_mover_interp(points1.select(startswith(emb))+1e-5,
                                   points2.select(startswith(emb))+1e-5,
                                   normer=normer,
                                  )**2
        #sums[0].append(points1.select(startswith(emb)).mean())
        #sums[1].append(points2.select(startswith(emb)).mean())
    #dist += earth_mover(np.array(sums[0]), np.array(sums[1]))
    return (dist/len(embs))**.5

def mp_earth_mover(args):
    i, j = args
    return earth_mover(i, j)

def mp_earth_mover_multi(args):
    i, j = args
    return earth_mover_multi(i, j)

import progressbar as pb
def pdist(X, metric, p=2, w=None, V=None, VI=None):
    X = np.asarray(X, order='c')


    s = X.shape
    if len(s) != 2:
        raise ValueError('A 2-dimensional array must be passed.')

    m, n = s
    dm = np.zeros((m * (m - 1) / 2,), dtype=np.double)



    k = 0
    prog = pb.ProgressBar(widgets=['calculating distances', pb.Bar(),
                                   pb.Percentage(), pb.ETA()])
    for i in prog(range(0, m - 1)):
        for j in range(i + 1, m):
            dm[k] = metric(X[i], X[j])
            k = k + 1
    prog.finish()
    return dm

def mp_mapped(args):
    manager, X, i, j = args
    metric = manager.get_metric()
    return metric(X[i], X[j])

def mp_pdist(X, metric, p=2, w=None, V=None, VI=None):
    import multiprocessing
    from multiprocessing.managers import BaseManager
    X = np.asarray(X, order='c')


    s = X.shape
    if len(s) != 2:
        raise ValueError('A 2-dimensional array must be passed.')

    m, n = s
    dm = np.zeros((m * (m - 1) / 2,), dtype=np.double)

    pool = multiprocessing.Pool(10)
    func = globals()["mp_"+metric.__name__]

    k = 0
    prog = pb.ProgressBar(widgets=['calculating distances', pb.Bar(),
                                   pb.Percentage(), pb.ETA()])
    for i in prog(range(0, m - 1)):
        ks = np.arange(k, k + m - i - 1)
        inputs = [(X[i], X[j]) for j in range(i+1, m)]
        dm[ks] = pool.map(func, inputs)
        k  = ks[-1] + 1
    prog.finish()
    pool.close()
    return dm

def mp_pandas_pdist(X, metric, p=2, w=None, V=None, VI=None):
    import multiprocessing

    s = X.shape
    if len(s) != 2:
        raise ValueError('A 2-dimensional array must be passed.')

    m, n = s
    dm = np.zeros((m * (m - 1) / 2,), dtype=np.double)

    pool = multiprocessing.Pool()
    if metric.__name__.endswith('multi'):
        func = globals()["mp_"+metric.__name__]
    else:
        func = globals()["mp_"+metric.__name__+"_multi"]

    k = 0
    prog = pb.ProgressBar(widgets=['calculating distances', pb.Bar(),
                                   pb.Percentage(), pb.ETA()])
    for i in prog(range(0, m - 1)):
        ks = np.arange(k, k + m - i - 1)
        inputs = [(X.ix[i], X.ix[j]) for j in range(i+1, m)]
        dm[ks] = pool.map(func, inputs)
        k  = ks[-1] + 1
    prog.finish()
    pool.close()
    return dm

def pandas_pdist(X, metric, p=2, w=None, V=None, VI=None):
    s = X.shape
    if len(s) != 2:
        raise ValueError('A 2-dimensional array must be passed.')

    m, n = s
    dm = np.zeros((m * (m - 1) / 2,), dtype=np.double)

    if metric.__name__.endswith('multi'):
        func = globals()["mp_"+metric.__name__]
    else:
        func = globals()["mp_"+metric.__name__+"_multi"]

    k = 0
    prog = pb.ProgressBar(widgets=['calculating distances', pb.Bar(),
                                   pb.Percentage(), pb.ETA()])
    for i in prog(range(0, m - 1)):
        ks = np.arange(k, k + m - i - 1)
        inputs = [(X.ix[i], X.ix[j]) for j in range(i+1, m)]
        print(len(inputs))
        print(ks)
        dm[ks] = list(map(func, inputs))
        k  = ks[-1] + 1
    prog.finish()
    return dm


if __name__ == "__main__":
    import pandas as pd
    import matplotlib.pyplot as mpl
    from multiprocessing import Pool

    eps = 1.01
    kwargs = dict(index_col=0,
            keep_default_na=False,
            na_values=['---'])
    ase = pd.read_table('analysis_godot/ase_summary.tsv', **kwargs)
    expr = pd.read_table('analysis_godot/summary_fb.tsv', **kwargs)
    translate = pd.read_table('prereqs/gene_map_table_fb_2016_01.tsv', index_col=1).ix[:,0]
    with Pool() as p:
        diff_from_uniform = pd.Series(
                data=map(
                    get_nonuniform_mp,
                    ase.iterrows()),
                index=ase.index
                )
        diff_from_uniform2 = pd.Series(
                data=map(
                    get_nonuniform_inv_mp,
                    ase.iterrows()),
                index=ase.index
                )

    good_ase = pd.np.isfinite(ase).sum(axis=1) > 25
    lott = pd.read_table('prereqs/journal.pbio.1000590.s002', index_col=0)
    is_mat = {gene for gene, c in lott.CLASS.iteritems() if c == 'mat'}
    is_mat = translate.apply(is_mat.__contains__)

    diff_from_uniform = diff_from_uniform.ix[~(is_mat.ix[diff_from_uniform.index].replace(pd.np.nan, False)) & good_ase].sort_values(ascending=False)
    diff_from_uniform2 = diff_from_uniform2.ix[~(is_mat.ix[diff_from_uniform2.index].replace(pd.np.nan, False)) & good_ase].sort_values(ascending=False)
    diff_from_uniform.to_csv('analysis/results/diff_from_uniform.tsv', sep='\t')
    diff_from_uniform2.to_csv('analysis/results/diff_from_uniform2.tsv', sep='\t')


    import PlotUtils


    ix = diff_from_uniform.index[:50].intersection(expr.index)
    ix2 = diff_from_uniform2.index[:50].intersection(expr.index)
    plot_kwargs = dict(
            draw_row_labels=True,
            box_size=15,total_width=150,
            split_columns=True, col_sep='_sl',
            convert=True,
            progress_bar=True,
            )

    PlotUtils.svg_heatmap(
            ase.ix[ix, :-1],
            'analysis/results/diff_from_uniform.svg',
            row_labels=translate.ix[ix],
            norm_rows_by='center0pre',
            cmap=mpl.cm.RdBu,
            **plot_kwargs
            )
    PlotUtils.svg_heatmap(
            expr.ix[ix,:-1],
            'analysis/results/diff_from_uniform_expr.svg',
            row_labels=translate.ix[ix],
            norm_rows_by='max',
            cmap=PlotUtils.ISH,
            **plot_kwargs
            )
    PlotUtils.svg_heatmap(
            ase.ix[ix2, :-1],
            'analysis/results/diff_from_uniform2.svg',
            row_labels=translate.ix[ix2],
            norm_rows_by='center0pre',
            cmap=mpl.cm.RdBu,
            **plot_kwargs
            )
    PlotUtils.svg_heatmap(
            expr.ix[ix2,:-1],
            'analysis/results/diff_from_uniform2_expr.svg',
            row_labels=translate.ix[ix2],
            norm_rows_by='max',
            cmap=PlotUtils.ISH,
            **plot_kwargs
            )

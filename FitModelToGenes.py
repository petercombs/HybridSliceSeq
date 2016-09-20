from __future__ import print_function
import PointClouds as pc
from sklearn.linear_model import logistic
from statsmodels.api import Logit
import PlotUtils as pu
from matplotlib import pyplot as mpl
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import numpy as np
import pandas as pd
import multiprocessing as mp
import itertools as it
from progressbar import ProgressBar as pbar
from collections import Counter


target_gene = 'hb'
input_genes = ['bcdP', 'kni', 'gtP']

all_regs = set(('brk bun cad CG10924 CG17786 CG4702 cnc croc hkb kni tll Traf1 trn '
            'tsh twi bcdP gtP bcdP2').split())

all_regs = set(('bcdP bcdP2 cad gtP kni KrP hkb hkb2 tll D da dl mad med shn '
                'sna twi').split())
all_regs = set(('bcdP bcdP2 cad gtP KrP tll D da dl mad med shn '
                'sna twi').split())


time_point = 'T2'

all_changes = {
    'hb': {
        'hkb' : 0.5, 'hkb2': 0.5, 'kni': 0.5, 'D': 1.0, 'KrP': 2.0,
        'gtP': 0.5, 'tll': 1.0, 'bcdP2': 2.0, 'bcdP': 2,
    },
    'Kr': {
        'gtP': 0.5, 'tll': 2.0, 'hkb': 0.5, 'kni': 2.0, 'D': 0.5,
        'hbP': 2.0,

    },
    'Bsg25A': {
        'gtP': 2.0, 'hkb': 2.0, 'KrP': 2.0, 'hbP': 0.5,

    }
}


changes = all_changes[target_gene]

def plot_changes(model, changes, atlas, new_figure=True, xlims=None, ylims=None,
                hide_low_expr=False):
    X = atlas.ix[:, model.params.index]
    p1 = model.predict(X)
    fitter3 = Predictor(model, changes)
    p2 = fitter3.predict(X)
    if new_figure:
        mpl.figure()
    mpl.scatter(
        atlas.x,
        atlas.z,
        s=60,
        c=(p2 - p1) / (p2 + p1),
        vmin=-1, vmax=1,
        edgecolors=(0,0,0,0),
        cmap=pu.cm.RdBu,
    )
    if xlims is not None:
        mpl.xlim(*xlims)
    if ylims is not None:
        mpl.ylim(*ylims)
    mpl.xticks([])
    mpl.yticks([])
    if hide_low_expr is not False:
        mpl.scatter(
            atlas.x[(p1 < hide_low_expr) & (p2 < hide_low_expr)],
            atlas.z[(p1 < hide_low_expr) & (p2 < hide_low_expr)],
            s=60,
            c=(.5, .5, .5),
            edgecolors=(0, 0, 0, 0),
        )

    mpl.ylabel('predicted ASE data')
    mpl.title(str(changes))
    mpl.gca().set_aspect(1)

class Predictor():
    def __init__(self, model, changes):
        self.params = model.params.copy()
        for tf in changes:
            if tf in self.params:
                self.params[tf] *= changes[tf]

    def predict(self, X):
        return 1/(1 + np.exp(
            -np.dot(self.params, X.T)
        )
        )


def fit_model(X, y, co=0.1):
    sm = Logit((y.clip(0, 1)>co).astype(float), X.clip(0, 1), missing='drop')
    return sm.fit(disp=False)

if __name__ == "__main__":
    if 'pcr' not in locals():
        pcr = pc.PointCloudReader(open('prereqs/D_mel_wt__atlas_r2.vpc'))
        atlas_expr, atlas_coords = pcr.data_to_arrays()

    for gene in atlas_expr.major_axis:
        last_t = None
        for t in atlas_expr.minor_axis:
            if not atlas_expr.ix[:, gene, t].count() and last_t:
                atlas_expr.ix[:, gene, t] = atlas_expr.ix[:, gene, last_t]
            elif atlas_expr.ix[:, gene, t].count() and not last_t:
                last_t = t


    atlas_expr.ix[:, 'const', :] = 1
    atlas_expr.ix[:, 'bcdP2', :] = atlas_expr.ix[:, 'bcdP', :]**2
    atlas_expr.ix[:, 'hkb2', :] = atlas_expr.ix[:, 'hkb', :]**2
    X = atlas_expr.ix[:, input_genes, time_point].T
    #X['bcd2'] = X['bcdP']**2
    X['const'] = 1
    y = atlas_expr.ix[:, target_gene,time_point]
    co = .2

    x_coord = atlas_coords.ix[:, 'X', time_point]
    z_coord = atlas_coords.ix[:, 'Z', time_point]
    x_coord_scale = x_coord.copy()
    x_coord_scale = ((x_coord_scale - x_coord_scale.min()) /
                      (x_coord_scale.max() - x_coord_scale.min()))

    poly1 = Path([
        [-300, 100],  # Upper left
        [-75, 100], [-93, 0], [-138, -100],  # inter-stripe boundary
        [-300, -100], [-300, 100], # lower left to completion
    ])
    poly2 = Path([
        [-75, 100], [-93, 0], [-138, -100],  # inter-stripe boundary
        [58, -100], [28, 100], # inter-stripe boundary
        [-75, 100]
    ])

    poly12 = Path([
        [-230, 90],  # Upper left
        [-230, -100],
        [40, -100], [25, 90], # inter-stripe boundary
        [-230, 90],  # Upper left

    ])

    poly3 = Path([
        [58, -100], [28, 100], # inter-stripe boundary
        [300, 100], [300, -100], [58, -100],
    ])

    kr_poly = Path([
        [-68, 90], [70, 90],
        [91, -100], [-100, -100],
        [-68, 90],
    ])

    all_poly = Path([
        [-230, 90], [300, 90],
        [300, -100], [-230, -100],
        [-230, 90]
    ])

    small_atlas = pd.DataFrame({
            'x' : atlas_coords.ix[:, 'X', time_point],
            'z' : atlas_coords.ix[:, 'Z', time_point],
            'c' : atlas_expr.ix[:, target_gene, time_point],
    })

    polys = {
        'hb': [poly12],
        'Kr': [kr_poly],
        'Bsg25A': [all_poly],
    }

    for poly in polys[target_gene]:

        in_central = poly.contains_points(
            atlas_coords.ix[:, ['X', 'Z'], time_point].T
        )
        not_expr = atlas_expr.ix[:, target_gene, time_point] < co
        in_central |= not_expr
        print(sum(in_central))
        #in_central =  (x_coord < 45)
        #in_central = x_coord_scale < 0.6

        #fitter = logistic.LogisticRegression(fit_intercept=False)
        #fitter.fit(X.ix[in_central, :], y.ix[in_central] > co)

        sm_fitter = Logit( y.ix[in_central].clip(0, 1), X.ix[in_central].clip(0, 1))
        sm_fit = sm_fitter.fit()

        Y_tmp = atlas_expr.ix[in_central, target_gene,time_point].copy()

        all_regs = atlas_expr.ix[:, all_regs, time_point].count(axis=1) > 0
        all_regs = all_regs.index[all_regs]


        if True:
        #if (poly == poly1) or (poly == poly2) or (poly == poly12):
        #if False:
            #best_tfs = ['bcdP', 'hkb', 'hkb2', 'KrP', 'bcdP2', 'const']
            #best_tfs = ['bcdP', 'bcdP2', 'gtP', 'kni', 'hkb',  'KrP', 'const']
            #best_tfs = atlas_expr.major_axis
            best_tfs = {'bcdP', 'bcdP2', 'cad', 'gtP', 'hbP', 'KrP', 'hkb', 'tll', 'D',
                        'h', 'dl', 'mad', 'shn', 'twi', 'const'}
            best_tfs.discard(target_gene)
            best_tfs.discard(target_gene + 'P')
            best_model = fit_model(atlas_expr.ix[in_central, best_tfs,
                                                 time_point].dropna(how='all').T,
                                   Y_tmp, co=co)
            best_tfs = best_model.params.index
        else:
            with mp.Pool(1) as p:
                pool = {}
                for comb in it.combinations(all_regs, 4):
                    comb = comb + ('hkb', 'kni', 'const', )
                    X_tmp = atlas_expr.ix[in_central, comb, time_point].T.copy()
                    pool[comb] = p.apply_async(fit_model, (X_tmp, Y_tmp, co))
                outs = {}
                pr2 = pd.Series(index=pool.keys(), data=np.nan)
                llrs =pd.Series(index=pool.keys(), data=np.nan)
                for comb in pbar()(pool):
                    outs[comb] = pool[comb].get()
                    pr2[comb] = outs[comb].prsquared
                    llrs[comb] = outs[comb].llr
            best_tfs = pr2.sort_values().index[-1]
            best_model = outs[best_tfs]

        print(best_model.summary().as_text())
        best_X = atlas_expr.ix[:, best_tfs, time_point].T

        small_atlas['in_central'] = in_central
        small_atlas['color'] = [
            'b' if not ic else 'k' if yy > co else 'w'
            for ic, yy in zip(small_atlas.in_central, small_atlas.c)
        ]
        for tf in best_tfs:
            small_atlas[tf] = atlas_expr.ix[:, tf, time_point]
        small_atlas.sort_values(by='c', inplace=True)

        ax = mpl.subplot(3, 2, 1, aspect='equal')
        mpl.scatter(
            small_atlas.x,
            small_atlas.z,
            s=60,
            c=small_atlas.c,
            edgecolors=(0, 0, 0, 0),
            cmap=pu.ISH,
        )
        mpl.ylabel('Original data')
        print(len(small_atlas.c))
        ax.add_patch(PathPatch(poly, fill=False))
        xlims = mpl.xlim()
        xlims = (xlims[0] - 10, xlims[1] + 10)
        mpl.xlim(*xlims)
        ylims = mpl.ylim()
        ylims = (ylims[0] - 20, ylims[1] + 20)
        mpl.ylim(*ylims)
        mpl.yticks([])

        '''
        mpl.subplot(3, 2, 2, aspect='equal')
        mpl.scatter(
            small_atlas.x,
            small_atlas.z,
            s=60,
            c=small_atlas.c,
            edgecolors=small_atlas.color,
            cmap=pu.ISH,

        )
        mpl.yticks([])
        mpl.ylabel('Fitting')
        '''

        mpl.subplot(3, 2, 3, aspect='equal')

        p1 = best_model.predict(best_X.ix[in_central])
        mpl.scatter(
            atlas_coords.ix[in_central, 'X', time_point],
            atlas_coords.ix[in_central, 'Z', time_point],
            s=60,
            c=p1,
            edgecolors=(0,0,0,0),
            cmap=pu.ISH,
        )
        #print(fitter.score(X.ix[in_central, :], y.ix[in_central] > co))
        mpl.xlim(*xlims)
        mpl.yticks([])
        mpl.ylabel('Predicted Mel Data')

        mpl.subplot(3, 2, 5, aspect='equal')

        fitter2 = Predictor(best_model, changes)
        p2 = fitter2.predict(best_X.ix[in_central])
        mpl.scatter(
            atlas_coords.ix[in_central, 'X', time_point],
            atlas_coords.ix[in_central, 'Z', time_point],
            s=60,
            c=(p2 - p1)/(p2 + p1),
            vmin=-1, vmax=1,
            edgecolors=(0,0,0,0),
            cmap=pu.cm.RdBu,
        )
        mpl.xlim(*xlims)
        mpl.yticks([])
        mpl.ylabel('predicted ASE data')

        mpl.subplot(3, 2, 4, aspect='equal')
        mpl.scatter(
            atlas_coords.ix[in_central, 'X', time_point],
            atlas_coords.ix[in_central, 'Z', time_point],
            s=60,
            c=p2,
            edgecolors=(0,0,0,0),
            cmap=pu.ISH,
        )
        mpl.xlim(*xlims)
        mpl.yticks([])
        mpl.ylabel('predicted sim data')

    mpl.savefig('analysis/results/model_tweak/{}'.format(target_gene))
    print(best_model.summary().as_latex(),
          file=open('analysis/results/model_tweak/{}_model.tex'.format(target_gene), 'w')
         )
    for tf in (best_model
               .params.index.intersection(
                   [key for key, value in changes.items() if value != 1])
              ):
        plot_changes(best_model, {tf: changes[tf]}, small_atlas.ix[small_atlas.in_central], xlims=xlims,
                     ylims=ylims)
        mpl.savefig('analysis/results/model_tweak/{}_{}'
                    .format(target_gene, tf))

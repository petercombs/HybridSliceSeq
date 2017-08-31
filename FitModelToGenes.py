from __future__ import print_function
import PointClouds as pc
#from sklearn.linear_model import logistic
from statsmodels.api import Logit
import PlotUtils as pu
import matplotlib as mplbase
from matplotlib import pyplot as mpl
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.transforms import Bbox
import numpy as np
import pandas as pd
import multiprocessing as mp
import itertools as it
import setcolor
from progressbar import ProgressBar as pbar
from collections import defaultdict
from sys import argv
from os import path
import Utils as ut
from argparse import ArgumentParser
from re import sub


#target_gene = 'hb'
target_gene = 'comm2'
input_genes = ['bcdP', 'kni', 'gtP']

all_regs = set(('brk bun cad CG10924 CG17786 CG4702 cnc croc hkb kni tll Traf1 trn '
            'tsh twi bcdP gtP bcdP2').split())

all_regs = set(('bcdP bcdP2 cad gtP kni KrP hkb hkb2 tll D da dl mad med shn '
                'sna twi').split())
all_regs = set(('bcdP bcdP2 cad gtP KrP tll D da dl mad med shn '
                'sna twi zen brk emc numb rho tkv Doc2').split())


times = {
    'Ance': 'T4',
    'hb': 'T3',
    'path': 'T3',
    'Bsg25A': 'T4',
    'pxb': 'T4',
    'bmm': 'T4',
    'comm2': 'T4',
    'slp1': 'T2',
    'CG8147': 'T4',
    'sala': 'T4',
}

all_changes = {
    'Ance': {
        'gtP': 0.5, 'Kr': 2.0, 'zen': 2.0,
        'Doc2': 0.5, 'twi': 2.0, 'kni': 0.5,
        'tll': 2.0, 'hkb': 0.5, 'D': 0.5, 'sna': 2.0,
    },
    'hb': {
        'D': 2.0, 'twi': 2.0,
        'gtP': 0.5, 'tll': 1.0,
        'bcdP2': 2.0, 'bcdP': 2,
        'hkb' : 0.5, #'hkb2': 0.5,
        'kni': 0.5, 'D': 1.0, 'KrP': 0.5,
        'sna' : 2.0,
    },
    'Kr': {
        'gtP': 0.5, 'tll': 2.0, 'hkb': 0.5, 'kni': 2.0, 'D': 0.5,
        'hbP': 2.0,

    },
    'Bsg25A': {
        'gtP': 2.0, 'hkb': 2.0, 'KrP': 2.0, 'hbP': 0.5,

    },
    'path': {
        'twi': 0.5, 'hb': 0.5, 'kni': 2.0,
        'hkb': 2.0, 'cad': 0.5, 'tll': 2.0,

    },
    'comm2': {
        'hkb': 2.0, 'KrP': 0.5, 'gtP': 2.0,
        'D': 0.5,

    },
    'bmm':{
        'D': 0.5, 'cad': .5,
        'sna': 0.5,
        'Kr': 2.0, 'hb': 2.0, 'kni':2.0,
    },
    'CG8147':{
        'twi': 0.5, 'hkb': 2.0, 'Kr':2.0, 'D': 2.0,
        'tll': 0.5,
    },
    'pxb': {
        'D': 0.5, 'hb': 0.5, 'twi':0.5,
        'Kr': 2.0,
        'kni': 1.0, 'tll': 1.0,
    },
    'sala': {
        'Doc2': 1.5, 'Kr': 0.5, 'hkb': 1.0, 'cad': 1.0, 'kni': 1.0,
        'D': 2.0,

    }
}



def plot_changes(model, changes, atlas, new_figure=True, xlims=None, ylims=None,
                 simplify_title=False,
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

    #mpl.ylabel('predicted ASE data')
    if simplify_title:
        mpl.title(', '.join(['Increase ', 'Decrease '][value < 1] + key
                            for key, value in changes.items()))
    else:
        mpl.title(str(changes))
    ax = mpl.gca()
    ax.set_aspect(1)
    ax.set_axis_bgcolor((.6, .6, .6))
    return p1, p2

def get_virtual_ase_slices(model, changes, atlas, n_slices, fitter=None,
                           denom_cutoff=5):
    ''' Calculate predicted ASE with model and single TF changes

    To use a fixed value for expression levels in the denominator, use a
    positive value for denom_cutoff.

    To have a fixed number of slices with non-nan expression, use a negative
    value for denom_cutoff.

    '''
    X = atlas.ix[:, model.params.index]
    p1 = model.predict(X)
    if fitter is None:
        fitter = Predictor(model, changes)
    p2 = fitter.predict(X)
    virtualslices = np.empty((1, n_slices)) * np.nan
    denoms = np.zeros(n_slices)
    pct_x = atlas.x - atlas.x.min()
    pct_x *= 100 /pct_x.max()
    for i in range(n_slices):
        in_slice = np.array(((100.0 * i / n_slices < pct_x)
                    & (pct_x < 100.0 * (i+1) / n_slices)))
        mel_in_slice = p1[in_slice].clip(0, 1).sum()
        sim_in_slice = p2[in_slice].clip(0, 1).sum()
        denom_i = mel_in_slice + sim_in_slice
        if denom_i > denom_cutoff:
            virtualslices[0, i] = (sim_in_slice - mel_in_slice) / np.sqrt(denom_i)
        else:
            virtualslices[0, i] = np.nan
        denoms[i] = denom_i
    if denom_cutoff < 0:
        virtualslices[0, denoms < sorted(denoms)[denom_cutoff]] = np.nan
    return denoms, virtualslices, fitter


def compare_ase(ase, model, changes, atlas, denom_cutoff=5):
    adaptive = denom_cutoff == 'adaptive'
    embryos = {col.split('sl')[0] for col in ase.index}
    pred_ase = pd.Series(index=ase.index, data=np.nan)
    pred_expr = pd.Series(index=ase.index, data=np.nan)
    fitter = None
    for embryo in embryos:
        emb_ase = ase.select(ut.startswith(embryo))
        if adaptive:
            denom_cutoff = -emb_ase.count()
        (denoms, virtualase, fitter) = get_virtual_ase_slices(
            model, changes, atlas, len(emb_ase),
            fitter=fitter, denom_cutoff=denom_cutoff,
        )
        pred_ase.ix[emb_ase.index] = virtualase[0]
        pred_expr.ix[emb_ase.index] = denoms

    return pred_ase.corr(ase), pred_ase, pred_expr

def calculate_r2(ase, prediction, mean_norm = False):
    if mean_norm:
        ase = ase - ase.mean()
        prediction = prediction - prediction.mean()
    return 1-(
                    ((ase - prediction)**2).sum()
                    / ((ase - ase.mean())**2).sum()
                )

def compare_ase_range(ase, model, tfs, atlas, range_low=0.1, range_high=10,
                      num_points=51, spacing='log',
                      denom_cutoff=5):
    if spacing == 'log':
        values = np.exp(np.linspace(np.log(range_low), np.log(range_high), num_points,
                        endpoint=True))
    elif spacing == 'linear':
        values = np.linspace(range_low, range_high, num_points, endpoint=True)
    else:
        raise ValueError("argument 'spacing' should be one of {'log', 'linear'}")

    if not (isinstance(tfs, tuple) or isinstance(tfs, list)):
        tfs = [tfs]


    corrs = np.zeros_like(values)
    r2s = np.zeros_like(values)
    for i, value in enumerate(values):
        corr, pred_ase, pred_expr = compare_ase(ase, model, {tf: value for tf in tfs},
                                 atlas, denom_cutoff)
        corrs[i] = corr
        r2s[i] = calculate_r2(ase, pred_ase, mean_norm=True)
    return values, corrs, r2s




def plot_virtual_ase_slices(model, changes, atlas, new_figure=True, n_slices=25,
                            xlims=None,
                            fitter=None, denom_cutoff=5):

    denoms, virtualslices, fitter = get_virtual_ase_slices(
        model, changes, atlas, n_slices,
        fitter=fitter, denom_cutoff=denom_cutoff
    )
    if new_figure:
        mpl.figure()
    virtualslices = np.ma.masked_where(np.isnan(virtualslices), virtualslices)
    mpl.cm.RdBu.set_bad((0.5, 0.5, 0.5))
    mpl.cm.RdBu.set_over((0.5, 0.5, 0.5))
    mpl.cm.RdBu.set_under((0.5, 0.5, 0.5))
    xs = np.linspace(atlas.x.min(), atlas.x.max(), n_slices, endpoint=True)
    mpl.pcolormesh(xs, [0,1], virtualslices, cmap=mpl.cm.RdBu, vmin=-1, vmax=1)
    ax = mpl.gca()
    #ax.set_aspect(1)
    ax.hlines([0, 1], xs[0], xs[-1])
    ax.vlines([xs[0], xs[-1]], 0, 1)
    if xlims:
        ax.set_xlim(xlims)
    ax.set_ylim(-0.1, 1.1)
    ax.set_xticks([])
    ax.set_yticks([])
    return denoms, virtualslices

def plot_both(model, changes, atlas, n_slices=27, fitter=None, denom_cutoff=5,
              plot_other=2,
              xlims=None, ylims=None):
    fig = mpl.figure()
    width = .9
    left = 0.5 - width / 2
    width2 = 420/800 if plot_other else 650/800
    left2 = 0.5-width2 / 2
    height = 1/(2 + 1 + plot_other)-.01
    ax1 = fig.add_axes([left,height, width, 2*height], aspect=1,xlim=xlims)
    #ax1 = mpl.subplot2grid((3+plot_other,1), (plot_other,0), 2, 1)
    p1, p2 = plot_changes(model, changes, atlas.ix[atlas.in_central], xlims=xlims,
                          simplify_title=True,
                          ylims=ylims, new_figure=False)
    #ax2 = mpl.subplot2grid((3+plot_other,1), (2+plot_other,0), 1, 1)
    ax2 = fig.add_axes([left2, 0, width2, height])
    plot_virtual_ase_slices(model, changes, atlas.ix[atlas.in_central],
                            fitter=fitter,
                            xlims=xlims,
                            new_figure=False, n_slices=n_slices, denom_cutoff=denom_cutoff)
    if plot_other:
        #ax0 = mpl.subplot2grid((5, 1), (0, 0), plot_other, 1)
        ax0 = fig.add_axes([left, 3*height+.04, width, plot_other * height])
        ax0.scatter(atlas.x, atlas.z, s=60,
                    c=p2,
                    edgecolors=(0, 0, 0, 0),
                    cmap=pu.ISH,
                   )
        ax0.set_aspect(1)
        ax0.set_axis_bgcolor((.6, .6, .6))
        ax0.set_xticks([])
        ax0.set_yticks([])
        if xlims:
            ax0.set_xlim(*xlims)
        if ylims:
            ax0.set_ylim(*ylims)

    ax2.set_frame_on(False)
    #mpl.tight_layout()
    #mpl.ioff()
    #mpl.draw()
    #mpl.show(block=False)
    #mpl.ion()
    #mpl.draw_if_interactive()
    #mpl.show()
    #bbox1_pts = ax1.get_position().get_points()
    #ax1.set_position(Bbox(bbox1_pts))
    #mpl.draw()
    #mpl.show(block=False)
    #bbox1_pts = ax1.get_position().get_points()
    #ax1.set_position(Bbox(bbox1_pts))
    #bbox2_pts = ax2.get_position().get_points()
    #bbox2_pts[:, 0] = bbox1_pts[:, 0]
    #ax1.set_position(Bbox(bbox1_pts))
    #ax2.set_position(Bbox(bbox2_pts))
    return fig, ax1, ax2

def set_widths(ax1, ax2):
    bbox1_pts = ax1.get_position().get_points()
    bbox2_pts = ax2.get_position().get_points()
    bbox2_pts[:, 0] = bbox1_pts[:, 0]
    ax1.set_position(Bbox(bbox1_pts))
    ax2.set_position(Bbox(bbox2_pts))


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

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--atlas',
                        default='/home/pcombs/HybridSliceSeq/prereqs/D_mel_wt__atlas_r2.vpc')
    parser.add_argument('--dir', default='/home/pcombs/HybridSliceSeq/')
    parser.add_argument('--find-best-tweaks', default=False,
                        action='store_true')
    parser.add_argument('target_gene', default='hb')
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    target_gene = args.target_gene
    time_point = times[target_gene]
    changes = all_changes[target_gene]
    print('-'*32, target_gene, '-'*32, sep='\n')
    dir = args.dir
    if 'pcr' not in locals():
        pcr = pc.PointCloudReader(open(
            args.atlas
        ))
        atlas_expr, atlas_coords = pcr.data_to_arrays()
        atlas_expr.ix[:, 'const', :] = 1
        atlas_expr.ix[:, 'bcdP2', :] = atlas_expr.ix[:, 'bcdP', :]**2
        #atlas_expr.ix[:, 'hkb2', :] = atlas_expr.ix[:, 'hkb', :]**2

    for gene in atlas_expr.major_axis:
        last_t = None
        for t in atlas_expr.minor_axis:
            if not atlas_expr.ix[:, gene, t].count() and last_t:
                atlas_expr.ix[:, gene, t] = atlas_expr.ix[:, gene, last_t]
            elif atlas_expr.ix[:, gene, t].count() and not last_t:
                last_t = t


    X = atlas_expr.ix[:, input_genes, time_point].T
    #X['bcd2'] = X['bcdP']**2
    X['const'] = 1
    y = atlas_expr.ix[:, target_gene,time_point]
    co = .1

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

    polys = defaultdict(lambda : [all_poly])
    polys.update({
        'hb': [poly12],
        'Kr': [kr_poly],
        'Bsg25A': [all_poly],
        'Ance': [all_poly],
    })


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
        Y_tmp /= Y_tmp.max()
        Y_tmp = 1.0 * (Y_tmp > .5)

        all_regs = atlas_expr.ix[:, all_regs, time_point].count(axis=1) > 0
        all_regs = all_regs.index[all_regs]


        #if True:
        #if (poly == poly1) or (poly == poly2) or (poly == poly12):
        if target_gene == 'hb':
            #best_tfs = ['bcdP', 'hkb', 'hkb2', 'KrP', 'bcdP2', 'const']
            #best_tfs = ['bcdP', 'bcdP2', 'gtP', 'kni', 'hkb',  'KrP', 'const']
            #best_tfs = atlas_expr.major_axis
            best_tfs = {'bcdP', 'bcdP2', 'cad', 'gtP', 'hbP', 'KrP', 'hkb',
                        'tll', 'D',  'kni',
                        'sna',
                        'h', 'dl', 'mad', 'shn', 'twi', 'const'}
            #best_tfs = {'const', 'bcdP', 'bcdP2', 'D', 'twi'}
            best_tfs.discard(target_gene)
            best_tfs.discard(target_gene + 'P')
            best_model = fit_model(atlas_expr.ix[in_central, best_tfs,
                                                 time_point].dropna(how='all').T,
                                   Y_tmp, co=co)
            best_tfs = best_model.params.index
        else:
            with mp.Pool() as p:
                pool = {}
                for comb in it.combinations(all_regs, 4):
                    comb = tuple(set(comb) |
                                 set(all_changes[target_gene].keys()))
                    comb = comb + ('const', )
                    if 'bcdP' in comb and 'bcdP2' not in comb:
                        comb = comb + ('bcdP2', )

                    X_tmp = (atlas_expr.ix[in_central, comb, time_point]
                             .T.copy()
                             .dropna(how='all', axis=1)
                            )
                    comb = tuple(X_tmp.columns)
                    if comb in pool: continue
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

        ax = mpl.subplot(1, 2, 1, aspect='equal')
        mpl.scatter(
            small_atlas.x,
            small_atlas.z,
            s=60,
            c=small_atlas.c,
            edgecolors=(0, 0, 0, 0),
            cmap=pu.ISH,
        )
        ax.set_axis_bgcolor((.6, .6, .6))
        mpl.ylabel('Original data')
        #print(len(small_atlas.c))
        #ax.add_patch(PathPatch(poly, fill=False))
        #xlims = mpl.xlim()
        #xlims = (xlims[0] - 10, xlims[1] + 10)
        emb_width = small_atlas.x.max() - small_atlas.x.min()
        xlims = (small_atlas.x.min() - .1 * emb_width,
                 small_atlas.x.max() + .1 * emb_width)
        mpl.xlim(*xlims)
        ylims = mpl.ylim()
        ylims = (ylims[0] - 40, ylims[1] + 40)
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

        ax = mpl.subplot(1, 2, 2, aspect='equal')

        p1 = best_model.predict(best_X.ix[in_central])
        mpl.scatter(
            atlas_coords.ix[in_central, 'X', time_point],
            atlas_coords.ix[in_central, 'Z', time_point],
            s=60,
            c=p1,
            edgecolors=(0,0,0,0),
            cmap=pu.ISH,
            vmax=1.2,
        )
        ax.set_axis_bgcolor((.6, .6, .6))
        #print(fitter.score(X.ix[in_central, :], y.ix[in_central] > co))
        mpl.xlim(*xlims)
        mpl.ylim(*ylims)
        mpl.yticks([])
        mpl.ylabel('Predicted Mel Data')

        '''
        ax = mpl.subplot(3, 2, 5, aspect='equal')

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
        ax.set_axis_bgcolor((.6, .6, .6))
        mpl.xlim(*xlims)
        mpl.ylim(*ylims)
        mpl.yticks([])
        mpl.ylabel('predicted ASE data')

        ax = mpl.subplot(3, 2, 4, aspect='equal')
        mpl.scatter(
            atlas_coords.ix[in_central, 'X', time_point],
            atlas_coords.ix[in_central, 'Z', time_point],
            s=60,
            c=p2,
            edgecolors=(0,0,0,0),
            cmap=pu.ISH,
        )
        ax.set_axis_bgcolor((.6, .6, .6))
        mpl.xlim(*xlims)
        mpl.ylim(*ylims)
        mpl.yticks([])
        mpl.ylabel('predicted sim data')
        '''

    mpl.savefig(path.join(dir,
                          'analysis/results/model_tweak/{}'.format(target_gene)))
    setcolor.set_backgroundcolor(mpl.gca(), 'k')
    setcolor.set_foregroundcolor(mpl.gca(), 'w')
    mpl.savefig(path.join(dir,
                          'analysis/results/model_tweak/{}_bgblk'.format(target_gene)))
    print(best_model.summary().as_latex(),
          file=open(path.join(dir,
                              'analysis/results/model_tweak/{}_model.tex'.format(target_gene)),
                    'w')
         )
    ase = (
        pd.read_table('analysis_godot/ase_summary_by_read2.tsv', **ut.pd_kwargs)
        .ix[target_gene]
        .select(ut.startswith(('melXsim', 'simXmel')))
    )
    rnaseq_expr = (
        pd.read_table('analysis_godot/summary.tsv', **ut.pd_kwargs)
        .drop('---', axis=1, errors='ignore')
        .select(**ut.sel_startswith(('melXsim', 'simXmel')))
        .rename(columns=lambda x: x.replace('FPKM', 'ase_value'))
    )
    pu.svg_heatmap(ase,
                   'analysis/results/model_tweak/{}_ase.svg'.format(target_gene),
                   norm_rows_by='center0pre',
                   cmap=pu.cm.RdBu,
                   **pu.kwargs_single_column
                  )

    n_good = ase.select(ut.startswith('melXsim_cyc14C_rep1')).count()
    for tf in (best_model
               .params.index.intersection(
                   [key for key, value in changes.items() if value != 1])
              ):
        mpl.cm.RdBu.set_bad((0.5, 0.5, 0.5))
        mpl.cm.RdBu_r.set_bad((0.5, 0.5, 0.5))
        fig, ax1, ax2 = plot_both(best_model, {tf: changes[tf]},
                                  small_atlas.ix[small_atlas.in_central],
                                  xlims=xlims, ylims=ylims,
                                  denom_cutoff=0)


        fig.savefig(path.join(dir, 'analysis/results/model_tweak/{}_{}'
                    .format(target_gene, tf)))
        setcolor.set_backgroundcolor(ax1, 'k')
        setcolor.set_foregroundcolor(ax1, 'w')
        fig.savefig(path.join(dir, 'analysis/results/model_tweak/{}_{}_bgblk'
                    .format(target_gene, tf)))
        print(
            {tf: changes[tf]},
            compare_ase(ase, best_model, {tf: changes[tf]},
                        small_atlas.ix[small_atlas.in_central],
                        denom_cutoff='adaptive'
                       )[0]
             )
        mpl.close(fig)


    if args.find_best_tweaks:
        tweaked = {}
        for tf in pbar()(changes):
            deltas, corrs, r2s = compare_ase_range(
                ase, best_model,
                [tf],
                small_atlas.ix[small_atlas.in_central],
                denom_cutoff=0,
                num_points=51,
                range_low=.05,
                range_high=20,
            )
            if tf == 'bcdP2':
                tf = 'bcd$^2$'
            elif tf == 'bcdP':
                tf = 'bcd'
            tweaked[tf] = deltas, corrs, r2s

        if 'bcdP' in best_model.params and 'bcdP2' in best_model.params:
            deltas, corrs, r2s = compare_ase_range(
                ase, best_model,
                ['bcdP', 'bcdP2'],
                small_atlas.ix[small_atlas.in_central],
                denom_cutoff=0,
                num_points=51,
                range_low=.05,
                range_high=20,
            )
            tweaked['bcd + bcd$^2$'] = deltas, corrs, r2s

        corrfig = mpl.figure(figsize=(16, 6))
        ax1 = mpl.subplot(1, 2, 1)
        ax2 = mpl.subplot(1, 2, 2)
        kwargs =(
            mplbase.cycler('linestyle', ['-', '--', '-.', ':'])
            *mplbase.cycler('color', ['b', 'g', 'r', 'c', 'm', 'y', 'k'])
        )
        for (tf, (deltas, corrs, r2s)), kwarg in zip(sorted(tweaked.items()),
                                                     kwargs):
            ax1.semilogx(deltas, corrs, basex=2, label=tf, **kwarg)
            ax2.semilogx(deltas, r2s.clip(-0.1, 1)*100, basex=2, label=tf,
                         **kwarg)
            try:
                the_tf = sub('[\^\$ \+]', '', tf)
                best_corr = deltas[np.nanargmax(corrs)]
                best_r2 = deltas[np.nanargmax(r2s)]
                if tf == 'bcd':
                    tfchange = {'bcdP': best_corr}
                elif tf == 'bcd$^2$':
                    tfchange = {'bcdP2': best_corr}
                elif tf == 'bcd + bcd$^2$':
                    tfchange = {'bcdP': best_corr, 'bcdP2': best_corr}
                else:
                    tfchange = {tf: best_corr}

                fig, a1, a2 = plot_both(best_model, tfchange,
                                        small_atlas.ix[small_atlas.in_central],
                                        denom_cutoff = 0, xlims=xlims, ylims=ylims)
                print('{}/analysis/results/model_tweak/paramsearch/{}_corr_{}'
                        .format(dir, target_gene, the_tf)
                )
                fig.savefig('{}/analysis/results/model_tweak/paramsearch/{}_corr_{}'
                        .format(dir, target_gene, the_tf))
                #mpl.close(fig)
                for a in tfchange:
                    tfchange[a] = best_r2
                fig, a1, a2 = plot_both(best_model, tfchange,
                                        small_atlas.ix[small_atlas.in_central],
                                        denom_cutoff = 0, xlims=xlims, ylims=ylims)
                mpl.savefig('{}/analysis/results/model_tweak/paramsearch/{}_varexp_{}'
                        .format(dir, target_gene, the_tf))
                #mpl.close(fig)
                fig2 = mpl.figure()
                if the_tf in rnaseq_expr.index:
                    mpl.scatter(rnaseq_expr.ix[the_tf, ase.index], ase)
                    plotted = True
                elif the_tf[-1] in ('2', 'P') and the_tf[:-1] in rnaseq_expr.index:
                    mpl.scatter(rnaseq_expr.ix[the_tf[:-1], ase.index], ase)
                    plotted = True
                else:
                    plotted = False
                if plotted:
                    print("Saving", '{}/analysis/results/model_tweak/paramsearch/{}_tfase_{}'
                                .format(dir, target_gene, the_tf))

                    mpl.title('spearman {}, pearson {}'
                              .format(
                                  rnaseq_expr.ix[the_tf].corr(ase, method='spearman'),
                                  rnaseq_expr.ix[the_tf].corr(ase, method='pearson'),
                              )
                             )
                    mpl.xlabel('FPKM of {}'.format(the_tf))
                    mpl.ylabel('ASE')
                    mpl.savefig('{}/analysis/results/model_tweak/paramsearch/{}_tfase_{}'
                                .format(dir, target_gene, the_tf))
                    mpl.close(fig2)
                else:
                    print("Couldn't find ", the_tf)

            except AssertionError as err:
                print(err)
                pass
        ax2.legend(loc='lower left')
        ax1.vlines(1, -1, 1)
        ax2.vlines(1, *ax2.get_ybound())
        ax1.set_xlabel('TF Binding Strength Multiplier')
        ax2.set_xlabel('TF Binding Strength Multiplier')
        ax1.set_ylabel('Correlation with ASE')
        ax2.set_ylabel('Percent Variance Explained')
        corrfig.savefig('{}/analysis/results/model_tweak/paramsearch/{}_param_search.png'
                    .format(dir, target_gene))





    if 'bcdP' in changes and 'bcdP2' in changes:
        fig = mpl.figure()
        ax1 = mpl.subplot2grid((3,1), (0,0), 2, 1)
        plot_changes(best_model, {'bcdP': 2.0, 'bcdP2': 2.0}, small_atlas.ix[small_atlas.in_central], xlims=xlims,
                     ylims=ylims, new_figure=False)
        ax2 = mpl.subplot2grid((3,1), (2,0), 2, 1)
        plot_virtual_ase_slices(best_model, {'bcdP': 1.5, 'bcdP2': 2.0}, small_atlas.ix[small_atlas.in_central],
                     new_figure=False, n_slices=27, denom_cutoff=-n_good)
        #mpl.tight_layout()
        mpl.tight_layout()
        mpl.savefig(path.join(dir, 'analysis/results/model_tweak/{}_bcdbcd2'
                    .format(target_gene)))
        mpl.ioff()
        mpl.draw()
        #mpl.show()
        mpl.ion()
        ax2.set_frame_on(False)
        bbox1_pts = ax1.get_position().get_points()
        bbox2_pts = ax2.get_position().get_points()
        bbox2_pts[:, 0] = bbox1_pts[:, 0]
        bbox1_pts[:, 0] = bbox2_pts[:, 0]
        ax1.set_position(Bbox(bbox1_pts))
        ax2.set_position(Bbox(bbox2_pts))
        mpl.savefig(path.join(dir, 'analysis/results/model_tweak/{}_bcdbcd2'
                    .format(target_gene)))
        setcolor.set_backgroundcolor(ax1, 'k')
        setcolor.set_foregroundcolor(ax1, 'w')
        mpl.savefig(path.join(dir, 'analysis/results/model_tweak/{}_bcdbcd2_bgblk'
                    .format(target_gene)))

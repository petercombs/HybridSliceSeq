from __future__ import print_function
import PointClouds as pc
from matplotlib.pyplot import (scatter, gca, figure, pcolormesh)
from matplotlib import cm
from numpy import zeros, zeros_like, nanmedian, nanpercentile
import numpy as np
import pandas as pd
from sys import stderr
from multiprocessing import Pool
from progressbar import ProgressBar

def find_best_match(s1_pos, s1_expr, s2_pos, s2_expr):
    """Find the best-matching nucleus within a given time point

    Note that unlike in Fowlkes 2011, I am using only the provided time point
    (and assuming that the input atlases have already been sliced), since the
    positions change, and so I don't know what to make of the "30 nearest nuclei
    to be matched" statement.

    """
    best_N = 30
    dists = (
        (s1_pos.pctX - s2_pos.pctX)**2
        + (s1_pos.pctY - s2_pos.pctY)**2
        + (s1_pos.pctZ - s2_pos.pctZ)**2
    )

    dists.sort_values(inplace=True)
    dists = dists[:best_N]

    expr_dists = ((s1_expr - s2_expr.ix[dists.index, :])**2).sum(axis=1)


    return expr_dists.idxmin()

def find_all_matches(s1_pos, s1_expr, s2_pos, s2_expr, pool=False, drop_genes=False):
    if 'X' in s1_pos.index:
        s1_pos = s1_pos.T
    if 'X' in s2_pos.index:
        s2_pos = s2_pos.T

    if 'bcd' in s1_expr.index:
        s1_expr = s1_expr.T
    if 'bcd' in s2_expr.index:
        s2_expr = s2_expr.T

    in_both = s1_expr.columns.intersection(s2_expr.columns)
    if drop_genes:
        in_both = in_both.drop(drop_genes)
    s1_expr = s1_expr.ix[:, in_both]
    s2_expr = s2_expr.ix[:, in_both]

    out = pd.Series(index=s1_expr.index, data=0)
    pbar = ProgressBar(maxval=len(out))
    if pool is None:
       pool = Pool()
    if pool is False:
        for nuc in pbar(out.index):
            out.ix[nuc] = (
                find_best_match(s1_pos.ix[nuc, :], s1_expr.ix[nuc,:], s2_pos, s2_expr)
            )
        return out

    results = []
    for nuc in out.index:
        results.append(
            pool.apply_async(
                find_best_match,
                (s1_pos.ix[nuc, :], s1_expr.ix[nuc, :], s2_pos, s2_expr)
            )
        )
    for res, ix in pbar(zip(results, out.index)):
        out.ix[ix] = res.get()

    return out




if __name__ == "__main__":
    if 'sim_atlas' not in locals():
        sim_atlas = pc.PointCloudReader(open('prereqs/dsim-20120727-r2-ZW.vpc'))
        mel_atlas = pc.PointCloudReader(open('prereqs/D_mel_wt__atlas_r2.vpc'))
    else:
        print("Using preloaded data", file=stderr)
    sim_atlas_expr, sim_atlas_pos = sim_atlas.data_to_arrays()
    mel_atlas_expr, mel_atlas_pos = mel_atlas.data_to_arrays()

    sim_atlas_pos.ix[:, 'pctX', :] = 100*(sim_atlas_pos[:, 'X', :].subtract(sim_atlas_pos[:, 'X', :].min(axis=1), axis=0)).divide(sim_atlas_pos[:, 'X', :].max(axis=1) - sim_atlas_pos[:, 'X', :].min(axis=1), axis=0)
    sim_atlas_pos.ix[:, 'pctY', :] = 100*(sim_atlas_pos[:, 'Y', :].subtract(sim_atlas_pos[:, 'Y', :].min(axis=1), axis=0)).divide(sim_atlas_pos[:, 'Y', :].max(axis=1) - sim_atlas_pos[:, 'Y', :].min(axis=1), axis=0)
    sim_atlas_pos.ix[:, 'pctZ', :] = 100*(sim_atlas_pos[:, 'Z', :].subtract(sim_atlas_pos[:, 'Z', :].min(axis=1), axis=0)).divide(sim_atlas_pos[:, 'Z', :].max(axis=1) - sim_atlas_pos[:, 'Z', :].min(axis=1), axis=0)
    mel_atlas_pos.ix[:, 'pctX', :] = 100*(mel_atlas_pos[:, 'X', :].subtract(mel_atlas_pos[:, 'X', :].min(axis=1), axis=0)).divide(mel_atlas_pos[:, 'X', :].max(axis=1) - mel_atlas_pos[:, 'X', :].min(axis=1), axis=0)
    mel_atlas_pos.ix[:, 'pctY', :] = 100*(mel_atlas_pos[:, 'Y', :].subtract(mel_atlas_pos[:, 'Y', :].min(axis=1), axis=0)).divide(mel_atlas_pos[:, 'Y', :].max(axis=1) - mel_atlas_pos[:, 'Y', :].min(axis=1), axis=0)
    mel_atlas_pos.ix[:, 'pctZ', :] = 100*(mel_atlas_pos[:, 'Z', :].subtract(mel_atlas_pos[:, 'Z', :].min(axis=1), axis=0)).divide(mel_atlas_pos[:, 'Z', :].max(axis=1) - mel_atlas_pos[:, 'Z', :].min(axis=1), axis=0)


    mel_front_10 = mel_atlas_pos.ix[:, 'pctX', -2] <= 10
    sim_front_10 = sim_atlas_pos.ix[:, 'pctX', -2] <= 10

    mel_selector = (mel_atlas_expr[:, 'hb', -2] > .3) & (mel_atlas_pos.ix[:, 'pctX', -2] > 75)
    sim_selector = (sim_atlas_expr[:, 'hb', -2] > .3) & (sim_atlas_pos.ix[:, 'pctX', -2] > 75)

    mel_background = ((50 < mel_atlas_pos.ix[:, 'pctX', -2])
                      & (mel_atlas_pos.ix[:, 'pctX', -2] < 75))
    sim_background = ((50 < sim_atlas_pos.ix[:, 'pctX', -2])
                      & (sim_atlas_pos.ix[:, 'pctX', -2] < 75))

    print(
        mel_atlas_expr.ix[mel_selector, 'hb', -2].mean(),
        sim_atlas_expr.ix[sim_selector, 'hb', -2].mean(),
        mel_atlas_expr.ix[mel_front_10, 'hb', -2].mean(),
        sim_atlas_expr.ix[sim_front_10, 'hb', -2].mean(),
    )

    scatter(mel_atlas_pos[:, 'X', 4], mel_atlas_pos[:, 'Z', 4], c=mel_selector, s=80)
    scatter(sim_atlas_pos[:, 'X', 4], sim_atlas_pos[:, 'Z', 4], c=sim_selector, s=80)

    target = 'hb'
    factor = 3
    mel_grid_expr = zeros((101//factor+1, 101//factor+1))
    mel_grid_ns = zeros_like(mel_grid_expr)
    sim_grid_expr = zeros((101//factor+1, 101//factor+1))
    sim_grid_ns = zeros_like(sim_grid_expr)

    mel_expr_at_stage = (
        mel_atlas_expr.ix[:, target, -2]
        - np.nanmean(mel_atlas_expr.ix[mel_background, target, -2])
    )#.clip(0, np.inf)
    mel_expr_at_stage /= np.nanpercentile(
        mel_atlas_expr.ix[~mel_background, target, -2],
        90
    )
    sim_expr_at_stage = (
        sim_atlas_expr.ix[:, target, -2]
        - np.nanmean(sim_atlas_expr.ix[sim_background, target, -2])
    )#.clip(0, np.inf)
    sim_expr_at_stage /= np.nanpercentile(
        sim_atlas_expr.ix[~sim_background, target, -2],
        90
    )
    for x, z, expr in zip(mel_atlas_pos.ix[:, 'pctX', -2],
                          mel_atlas_pos.ix[:, 'pctZ', -2],
                          mel_expr_at_stage):
        x = int(x)//factor
        z = int(z)//factor
        mel_grid_expr[z, x] += expr
        mel_grid_ns[z, x] += 1

    for x, z, expr in zip(sim_atlas_pos.ix[:, 'pctX', -2],
                          sim_atlas_pos.ix[:, 'pctZ', -2],
                          sim_expr_at_stage):
        x = int(x)//factor
        z = int(z)//factor
        sim_grid_expr[z, x] += expr
        sim_grid_ns[z, x] += 1

    mel_grid_avg = mel_grid_expr / mel_grid_ns
    sim_grid_avg = sim_grid_expr / sim_grid_ns


    mel_grid_lo = nanpercentile(mel_grid_avg, 10)
    mel_grid_hi = nanpercentile(mel_grid_avg, 90)
    sim_grid_lo = nanpercentile(sim_grid_avg, 10)
    sim_grid_hi = nanpercentile(sim_grid_avg, 90)
    mel_grid_avg = (mel_grid_avg - mel_grid_lo) / (mel_grid_hi - mel_grid_lo)
    sim_grid_avg = (sim_grid_avg - sim_grid_lo) / (sim_grid_hi - sim_grid_lo)

    ax = gca()
    ax.set_aspect(1)

    denom = (sim_grid_avg + mel_grid_avg)


    figure()
    heatmap = pcolormesh((sim_grid_avg - mel_grid_avg) / np.where(denom > .5,
                                                                  denom,  np.nan),
               vmin=-1, vmax=1, cmap=cm.RdBu)
    heatmap.cmap.set_bad((.5, .5, .5))
    heatmap.cmap.set_under((.5, .5, .5))

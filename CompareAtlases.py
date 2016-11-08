from __future__ import print_function
import PointClouds as pc
from matplotlib.pyplot import (scatter, gca, figure, pcolormesh)
from matplotlib import cm
from numpy import zeros, zeros_like, nanmedian, nanpercentile
import numpy as np
from sys import stderr

if __name__ == "__main__":
    if 'sim_atlas' not in locals():
        sim_atlas = pc.PointCloudReader(open('prereqs/dsim-20120727-r2-ZW.vpc'))
        mel_atlas = pc.PointCloudReader(open('prereqs/D_mel_wt__atlas_r2.vpc'))
    else:
        print("Using preloaded data", file=stderr)
    sim_atlas_expr, sim_atlas_pos = sim_atlas.data_to_arrays()
    mel_atlas_expr, mel_atlas_pos = mel_atlas.data_to_arrays()

    sim_atlas_pos.ix[:, 'pctX', :] = 100*(sim_atlas_pos[:, 'X', :].subtract(sim_atlas_pos[:, 'X', :].min(axis=1), axis=0)).divide(sim_atlas_pos[:, 'X', :].max(axis=1) - sim_atlas_pos[:, 'X', :].min(axis=1), axis=0)
    sim_atlas_pos.ix[:, 'pctZ', :] = 100*(sim_atlas_pos[:, 'Z', :].subtract(sim_atlas_pos[:, 'Z', :].min(axis=1), axis=0)).divide(sim_atlas_pos[:, 'Z', :].max(axis=1) - sim_atlas_pos[:, 'Z', :].min(axis=1), axis=0)
    mel_atlas_pos.ix[:, 'pctX', :] = 100*(mel_atlas_pos[:, 'X', :].subtract(mel_atlas_pos[:, 'X', :].min(axis=1), axis=0)).divide(mel_atlas_pos[:, 'X', :].max(axis=1) - mel_atlas_pos[:, 'X', :].min(axis=1), axis=0)
    mel_atlas_pos.ix[:, 'pctZ', :] = 100*(mel_atlas_pos[:, 'Z', :].subtract(mel_atlas_pos[:, 'Z', :].min(axis=1), axis=0)).divide(mel_atlas_pos[:, 'Z', :].max(axis=1) - mel_atlas_pos[:, 'Z', :].min(axis=1), axis=0)


    mel_front_10 = mel_atlas_pos.ix[:, 'pctX', -2] <= 10
    sim_front_10 = sim_atlas_pos.ix[:, 'pctX', -2] <= 10

    mel_selector = (mel_atlas_expr[:, 'hb', -2] > .3) & (mel_atlas_pos.ix[:, 'pctX', -2] > 75)
    sim_selector = (sim_atlas_expr[:, 'hb', -2] > .3) & (sim_atlas_pos.ix[:, 'pctX', -2] > 75)

    print(
        mel_atlas_expr.ix[mel_selector, 'hb', -2].mean(),
        sim_atlas_expr.ix[sim_selector, 'hb', -2].mean(),
        mel_atlas_expr.ix[mel_front_10, 'hb', -2].mean(),
        sim_atlas_expr.ix[sim_front_10, 'hb', -2].mean(),
    )

    scatter(mel_atlas_pos[:, 'X', 4], mel_atlas_pos[:, 'Z', 4], c=mel_selector, s=80)
    scatter(sim_atlas_pos[:, 'X', 4], sim_atlas_pos[:, 'Z', 4], c=sim_selector, s=80)

    factor = 3
    mel_grid_expr = zeros((101//factor+1, 101//factor+1))
    mel_grid_ns = zeros_like(mel_grid_expr)
    sim_grid_expr = zeros((101//factor+1, 101//factor+1))
    sim_grid_ns = zeros_like(sim_grid_expr)

    for x, z, expr in zip(mel_atlas_pos.ix[:, 'pctX', -2],
                          mel_atlas_pos.ix[:, 'pctZ', -2],
                          mel_atlas_expr.ix[:, 'hb', -2]):
        x = int(x)//factor
        z = int(z)//factor
        mel_grid_expr[z, x] += expr
        mel_grid_ns[z, x] += 1

    for x, z, expr in zip(sim_atlas_pos.ix[:, 'pctX', -2],
                          sim_atlas_pos.ix[:, 'pctZ', -2],
                          sim_atlas_expr.ix[:, 'hb', -2]):
        x = int(x)//factor
        z = int(z)//factor
        sim_grid_expr[z, x] += expr
        sim_grid_ns[z, x] += 1

    mel_grid_avg = mel_grid_expr / mel_grid_ns
    sim_grid_avg = sim_grid_expr / sim_grid_ns


    mel_grid_avg /= nanpercentile(mel_grid_avg, 90)
    sim_grid_avg /= nanpercentile(sim_grid_avg, 90)

    ax = gca()
    ax.set_aspect(1)

    denom = (sim_grid_avg + mel_grid_avg)


    figure()
    heatmap = pcolormesh((sim_grid_avg - mel_grid_avg) / np.where(denom > 1,
                                                                  denom,  np.nan),
               vmin=-1, vmax=1, cmap=cm.RdBu)
    heatmap.cmap.set_bad((.5, .5, .5))
    heatmap.cmap.set_under((.5, .5, .5))

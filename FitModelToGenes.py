import PointClouds as pc
from sklearn.linear_model import logistic
import PlotUtils as pu
from matplotlib import pyplot as mpl
import numpy as np


target_gene = 'Kr'
input_genes = ['bcdP', 'kni', 'hbP', 'gtP']
time_point = 'T3'

if __name__ == "__main__":
    if 'pcr' not in locals():
        pcr = pc.PointCloudReader(open('prereqs/D_mel_wt__atlas_r2.vpc'))
        atlas_expr, atlas_coords = pcr.data_to_arrays()

    X = atlas_expr.ix[:, input_genes, time_point].T
    y = atlas_expr.ix[:, target_gene, time_point]
    co = .2

    x_coord = atlas_coords.ix[:, 'X', time_point]
    in_central =  (x_coord < 45)

    fitter = logistic.LogisticRegression()
    fitter.fit(X.ix[in_central, :], y.ix[in_central] > co)

    mpl.subplot(4, 1, 1)
    mpl.scatter(
        atlas_coords.ix[:, 'X', time_point],
        atlas_coords.ix[:, 'Z', time_point],
        s=60,
        c=atlas_expr.ix[:, target_gene, time_point],
        edgecolors=['k' if yy > co else 'w' for yy in y],
        cmap=pu.ISH,
           )
    xlims = mpl.xlim()
    mpl.yticks([])
    mpl.ylabel('Fitting')
    mpl.gca().set_aspect(1)

    mpl.subplot(4, 1, 2)
    p1 = fitter.predict_proba(X.ix[in_central, :])[:, 1]
    mpl.scatter(
        atlas_coords.ix[in_central, 'X', time_point],
        atlas_coords.ix[in_central, 'Z', time_point],
        s=60,
        c=p1,
        edgecolors=(0,0,0,0),
        cmap=pu.ISH,
    )
    print(fitter.score(X.ix[in_central, :], y.ix[in_central] > co))
    mpl.xlim(*xlims)
    mpl.yticks([])
    mpl.ylabel('Original Mel Data')
    mpl.gca().set_aspect(1)

    mpl.subplot(4, 1, 3)

    fitter2 = logistic.LogisticRegression()
    fitter2.coef_ = np.copy(fitter.coef_)
    fitter2.intercept_ = fitter.intercept_
    #fitter2.coef_[0,0] += 4.5
    fitter2.coef_[0,2] += 8.5
    fitter2.intercept_ -= 5.5
    p2 = fitter2.predict_proba(X.ix[in_central, :])[:, 1]
    mpl.scatter(
        atlas_coords.ix[in_central, 'X', time_point],
        atlas_coords.ix[in_central, 'Z', time_point],
        s=60,
        c=p2 - p1,
        vmin=-1, vmax=1,
        edgecolors=(0,0,0,0),
        cmap=pu.cm.RdBu,
    )
    mpl.xlim(*xlims)
    mpl.yticks([])
    mpl.ylabel('predicted ASE data')
    mpl.gca().set_aspect(1)

    mpl.subplot(4, 1, 4)
    mpl.scatter(
        atlas_coords.ix[in_central, 'X', time_point],
        atlas_coords.ix[in_central, 'Z', time_point],
        s=60,
        c=p2,
        edgecolors=(0,0,0,0),
        cmap=pu.ISH,
    )
    mpl.xlim(*xlims)
    mpl.gca().set_aspect(1)
    mpl.yticks([])
    mpl.ylabel('predicted sim data')

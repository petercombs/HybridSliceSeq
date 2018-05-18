import pandas as pd
import numpy as np
import Utils as ut
import matplotlib.pyplot as mpl

import warnings
warnings.filterwarnings('ignore', category=FutureWarning)

if __name__ == "__main__":
    if 'ase' not in locals() or ('reload_ase' in locals() and locals()['reload_ase']):
        print("Reloading data")
        ase = (pd.read_table('analysis_godot/ase_summary_by_read.tsv', **ut.pd_kwargs)
               .dropna(how='all', axis=1)
               .dropna(how='all', axis=0)
               .select(**ut.sel_startswith(('melXsim', 'simXmel')))
              )
        all_ase = ase.copy()
        expr = (pd.read_table('analysis_godot/summary.tsv', **ut.pd_kwargs)
                .drop('---', axis=1, errors='ignore')
               .dropna(how='all', axis=1)
               .dropna(how='all', axis=0)
                #.dropna(how='all', axis=1)
               )
        reload_ase = False

        syns = ut.get_synonyms()
        chrom_of = ut.get_chroms(syns)

        males = ('melXsim_cyc14C_rep3', 'simXmel_cyc14C_rep2')
        on_x = [chrom_of[gene] == 'X' for gene in ase.index]
        is_male = [col.startswith(males) for col in ase.columns]
        ase_nomaleX = ase.copy()
        ase_nomaleX.loc[on_x, is_male] = pd.np.nan
        ase = ase_nomaleX


    lott = pd.read_table('prereqs/journal.pbio.1000590.s002', index_col=0,
                         na_values='', keep_default_na=False)
    lott_zyg = ut.true_index(lott.CLASS == 'zyg')
    xs = (ut.get_xs(ase)+.01) * [(-1 if l.startswith('mel') else 1) for l in
                           ase.columns]
    sxs = xs.sort_values()
    sxs_ase = sxs
    ase_x_sorted = ase.loc[lott_zyg, sxs.index]
    ase_x_sorted = ase.loc[:, sxs.index]

    num_slices = len(sxs.index)
    figwidth, figheight = 6,9
    fig = mpl.figure(figsize=(figwidth,figheight))
    total_area = figwidth*figheight
    subplot_area = total_area / num_slices
    subplot_dim = np.sqrt(subplot_area)
    cols = figwidth // subplot_dim
    rows = int(np.ceil(num_slices / cols))
    print(rows, cols)

    offset = 0
    adj_ase_corrs = []
    for i, (col1, col2) in enumerate(zip(sxs.index, sxs.index[1:])):
        if col1.startswith('mel') and col2.startswith('sim'):
            offset = 1
            continue
        i -= offset
        ax = mpl.subplot(rows, cols, i+1)
        mpl.hist2d(ase_x_sorted[col1], ase_x_sorted[col2],
                   bins=np.arange(-1.05, 1.05, .05),
                   cmin=0, cmax=10,
                  )
        corr = ase_x_sorted[col1].corr(ase_x_sorted[col2])
        adj_ase_corrs.append(corr)
        mpl.text(-1, 1,
                "{:.02f}".format(corr),
                 fontdict={'size': 8, 'color': 'white'},
                horizontalalignment='left',
                verticalalignment='top')
        ax.set_aspect(1)
        #if i % cols != 0:
        ax.set_yticks([])
        #if i < (cols * (rows-1)):
        ax.set_xticks([])
    print("ASE Corrs:", np.mean(adj_ase_corrs), '+/-', np.std(adj_ase_corrs))
    fig.subplots_adjust(hspace=0.05, wspace=0.05,
                        left=0.025, right=0.975,
                        top=0.975, bottom=0.025)
    mpl.savefig('analysis/results/adj_slice_ase_corr', dpi=300)

    offset=0
    adj_ase_corrs = []
    fig = mpl.figure(figsize=(figwidth,figheight))
    ase_x_sorted = ase_x_sorted.ix[lott_zyg]
    for i, (col1, col2) in enumerate(zip(sxs.index, sxs.index[1:])):
        if col1.startswith('mel') and col2.startswith('sim'):
            offset = 1
            continue
        i -= offset
        ax = mpl.subplot(rows, cols, i+1)
        mpl.hist2d(ase_x_sorted[col1], ase_x_sorted[col2],
                   bins=np.arange(-1.05, 1.05, .05),
                   cmin=0, cmax=5,
                  )
        corr = ase_x_sorted[col1].corr(ase_x_sorted[col2])
        adj_ase_corrs.append(corr)
        mpl.text(-1, 1,
                "{:.02f}".format(corr),
                 fontdict={'size': 8, 'color': 'white'},
                horizontalalignment='left',
                verticalalignment='top')
        ax.set_aspect(1)
        #if i % cols != 0:
        ax.set_yticks([])
        #if i < (cols * (rows-1)):
        ax.set_xticks([])
    print("Zygotic ASE Corrs:", np.mean(adj_ase_corrs), '+/-', np.std(adj_ase_corrs))
    fig.subplots_adjust(hspace=0.05, wspace=0.05,
                        left=0.025, right=0.975,
                        top=0.975, bottom=0.025)
    mpl.savefig('analysis/results/adj_slice_ase_corr_zyg', dpi=300)




    fig = mpl.figure(figsize=(6,9))
    xs = (ut.get_xs(expr)+.01) * [(-1 if l.startswith('mel') else 1) for l in
                           expr.columns]
    sxs = xs.sort_values()
    expr_x_sorted = expr.loc[:, sxs.index]

    num_slices = len(sxs.index)
    figwidth, figheight = 6.5,9
    fig = mpl.figure(figsize=(figwidth,figheight))
    total_area = figwidth*figheight
    subplot_area = total_area / num_slices
    subplot_dim = np.sqrt(subplot_area)
    cols = figwidth // subplot_dim
    rows = int(np.ceil(num_slices / cols))

    adj_expr_corrs = []
    for i, (col1, col2) in enumerate(zip(sxs.index, sxs.index[1:])):
        ax = mpl.subplot(rows, cols, i+1)
        mpl.scatter(expr[col1]+1, expr[col2]+1, s=0.1)
        corr = expr[col1].add(1).apply(np.log10).corr(expr[col2].add(1).apply(np.log10))
        adj_expr_corrs.append(corr)
        mpl.text(2, 999,
                "{:.02f}".format(corr),
                horizontalalignment='left',
                verticalalignment='top')
        ax.set_aspect(1)
        ax.set_xlim(1, 1000)
        ax.set_ylim(1, 1000)
        ax.set_xscale('log', basex=10)
        ax.set_yscale('log', basey=10)
        if i % cols != 0:
            ax.set_yticks([])
        if i < (cols * (rows-1)):
            ax.set_xticks([])
    print("Expr Corrs: ", np.mean(adj_expr_corrs), '+/-', np.std(adj_expr_corrs))
    fig.subplots_adjust(hspace=0.05, wspace=0.05,
                        left=0.05, right=0.975,
                        top=0.975, bottom=0.025)
    mpl.savefig('analysis/results/adj_slice_expr_corr', dpi=300)


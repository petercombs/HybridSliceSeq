import pandas as pd
import Utils as ut
import matplotlib.pyplot as mpl

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
        to_gn = pd.read_table('prereqs/gene_map_table_fb_2016_01.tsv',
                index_col=1, skiprows=4).loc[:,0]
        to_fbgn = ut.get_synonyms()

        syns = ut.get_synonyms()
        chrom_of = ut.get_chroms(syns)

        males = ('melXsim_cyc14C_rep3', 'simXmel_cyc14C_rep2')
        on_x = [chrom_of[gene] == 'X' for gene in ase.index]
        is_male = [col.startswith(males) for col in ase.columns]
        ase_nomaleX = ase.copy()
        ase_nomaleX.loc[on_x, is_male] = pd.np.nan
        ase = ase_nomaleX


    xs = ut.get_xs(ase)
    sxs = xs.sort_values()
    ase_x_sorted = ase.loc[:, sxs.index]

    num_slices = len(sxs.index)
    figwidth, figheight = 6,9
    fig = figure(figsize=(figwidth,figheight))
    total_area = figwidth*figheight
    subplot_area = total_area / num_slices
    subplot_dim = np.sqrt(subplot_area)
    cols = figwidth // subplot_dim
    rows = int(np.ceil(num_slices / cols))

    for i, (col1, col2) in enumerate(zip(sxs.index, sxs.index[1:])):
        ax = mpl.subplot(rows, cols, i+1)
        mpl.scatter(ase[col1], ase[col2], s=0.1)
        mpl.text(-1, 1,
                "{:.02f}".format(ase[col1].corr(ase[col2])),
                horizontalalignment='left',
                verticalalignment='top')
    savefig('analysis/results/adj_slice_ase_corr', dpi=300)




    fig = figure(figsize=(6,9))
    xs = ut.get_xs(expr)
    sxs = xs.sort_values()
    expr_x_sorted = expr.loc[:, sxs.index]

import pandas as pd
import Utils as ut
import PlotUtils as pu
import numpy as np
from scipy import stats

if __name__ == "__main__":
    expr = locals()['expr']
    hyb = expr.select(**ut.sel_startswith(('melXsim', 'simXmel')))
    parental = expr.select(**ut.sel_startswith(('melXmel', 'simXsim')))
    melXmel = expr.select(**ut.sel_startswith('melXmel'))
    simXsim = expr.select(**ut.sel_startswith('simXsim'))

    is_mat = pd.read_table('analysis/results/maternal.tsv', squeeze=True,
                           header=None, **ut.pd_kwargs)
    mat_genes = is_mat.index[is_mat]

    xs = ut.get_xs(expr)
    hyb_xs = ut.get_xs(hyb)
    mel_xs = ut.get_xs(melXmel)
    sim_xs = ut.get_xs(simXsim)
    parental_xs = ut.get_xs(parental)

    kwargs = pu.kwargs_expr_heatmap.copy()
    kwargs['draw_row_labels'] = False
    kwargs['box_height'] = 1

    for tf, region in [('hb', (0, .1)),
                       ('Kr', (0.1, 0.4)),
                       ('Kr', (0.45, 0.75))]:
        peak_genes = {line.strip() for line in
                         open('analysis/results/{}_peak_genes_5k.txt'
                             .format(tf))
                         if line.strip() in expr.index}

        all_in_region = (region[0] < xs) & (xs < region[1])
        mel_in_region = (region[0] < mel_xs) & (mel_xs < region[1])
        sim_in_region = (region[0] < sim_xs) & (sim_xs < region[1])
        parental_in_region = (region[0] < parental_xs) & (parental_xs < region[1])

        gene_expr_level = parental.ix[peak_genes, parental_in_region].max(axis=1)
        expr_in_region = gene_expr_level.index[gene_expr_level > 5]

        non_mat_expr_genes = expr_in_region.difference(mat_genes)
        non_mat_expr_genes = gene_expr_level.ix[non_mat_expr_genes].sort_values().index

        pu.svg_heatmap(expr.ix[non_mat_expr_genes, all_in_region],
                       'analysis/results/non_mat_{}_{:0.2f}-{:0.2f}_expr_genes.svg'
                       .format(tf, region[0], region[1]),
                       squeeze_rows=np.nanmean,  progress_bar=True, **kwargs)


        mel_in_region = (melXmel.ix[non_mat_expr_genes, mel_in_region]
                         .divide(parental.ix[non_mat_expr_genes, parental_in_region]
                                 .max(axis=1), axis=0)
                         .mean(axis=1))
        sim_in_region = (simXsim.ix[non_mat_expr_genes, sim_in_region]
                         .divide(parental.ix[non_mat_expr_genes, parental_in_region]
                                 .max(axis=1), axis=0)
                         .mean(axis=1))
        res = stats.ttest_rel( mel_in_region, sim_in_region, nan_policy='omit',)

        print(tf, region, res, mel_in_region.mean(), sim_in_region.mean())


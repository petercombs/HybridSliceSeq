import pandas as pd
import Utils as ut
import PlotUtils as pu
import numpy as np
from scipy import stats
from GetASEStats import create_latex_command
from numpy import int8, int32, int64

if 'mpl' in locals():
    import matplotlib.pyplot as mpl
    HAS_MPL = True
else:
    HAS_MPL = False

EPSILON = .1

if __name__ == "__main__":
    data = {}
    expr = pd.read_table('godot/summary_wasp.tsv',
                         **ut.pd_kwargs).drop(['---'], axis=1)
    hyb = expr.select(**ut.sel_startswith(('melXsim', 'simXmel')))
    parental = expr.select(**ut.sel_startswith(('melXmel', 'simXsim')))
    melXmel = expr.select(**ut.sel_startswith('melXmel')).rename(columns=lambda x:x.split('_sl')[1])
    simXsim = expr.select(**ut.sel_startswith('simXsim')).rename(columns=lambda x:x.split('_sl')[1])
    higher_parental = melXmel.where(melXmel > simXsim, simXsim)
    lower_parental = melXmel.where(melXmel < simXsim, simXsim)
    mel_sim_ratio = (melXmel + EPSILON) / (simXsim + EPSILON)

    is_mat = pd.read_table('analysis/results/maternal.tsv', squeeze=True,
                           header=None, **ut.pd_kwargs)
    mat_genes = is_mat.index[is_mat]

    xs = ut.get_xs(expr)
    hyb_xs = ut.get_xs(hyb)
    mel_xs = np.linspace(0, 1, 27, endpoint=True)
    sim_xs = np.linspace(0, 1, 27, endpoint=True)
    parental_xs = ut.get_xs(parental)

    kwargs = pu.kwargs_expr_heatmap.copy()
    kwargs['draw_row_labels'] = False
    kwargs['box_height'] = 1
    kwargs['progress_bar'] = False

    distro_mels = []
    distro_sims = []

    target_regions = {
        'kr_ant_bg': ('Kr', (0.33, 0.44)),
        'Kr_central': ('Kr', (0.55, 0.72)),
        'Kr_post_bg': ('Kr', (0.74, 0.85)),
        'hb_ant_tip': ('hb', (0, .10)),
        'hb_ant_stripe': ('hb', (.20, .50)),
        'hb_interstripe': ('hb', (.60, .75)),
        'hb_post_stripe': ('hb', (.80, .90)),
        #('prd', (0, .1)),
    }

    nmegs = [] # Non maternal expressed gene lists

    for region_name, (tf, region) in target_regions.items():
        peak_genes = {line.strip() for line in
                         open(
                             #'analysis/results/{}_1000_peak_genes_2500bp_tss.txt'
                             'analysis/results/hb_wt_emd_0.1.txt'
                             #'analysis/results/{}_tss_genes.txt'
                             .format(tf.lower()))
                         if line.strip() in expr.index}
        data['{}diff'.format(tf)] = len(peak_genes)
        peak_genes.add(tf)

        all_in_region_ix = (region[0] <= xs) & (xs < region[1])
        mel_in_region_ix = (region[0] <= mel_xs) & (mel_xs < region[1])
        sim_in_region_ix = (region[0] <= sim_xs) & (sim_xs < region[1])
        parental_in_region_ix = (region[0] <= parental_xs) & (parental_xs < region[1])

        gene_expr_level = parental.ix[peak_genes, parental_in_region_ix].min(axis=1)
        expr_in_region = ut.true_index(gene_expr_level > -1)

        non_mat_expr_genes = expr_in_region.difference(mat_genes)
        non_mat_expr_genes = mel_sim_ratio.ix[non_mat_expr_genes,
                                              mel_in_region_ix].mean(axis=1).sort_values().index
        nmegs.append(non_mat_expr_genes)

        pu.svg_heatmap(expr.ix[non_mat_expr_genes, all_in_region_ix],
                       'analysis/results/non_mat_{}_{:0.2f}-{:0.2f}_expr_genes.svg'
                       .format(tf, region[0], region[1]),
                       squeeze_rows=np.nanmean,  **kwargs)


        mel_in_region = (melXmel.ix[non_mat_expr_genes, mel_in_region_ix]
                         .divide(parental.ix[non_mat_expr_genes, :]
                                 .max(axis=1), axis=0)
                         .mean(axis=1))
        sim_in_region = (simXsim.ix[non_mat_expr_genes, sim_in_region_ix]
                         .divide(parental.ix[non_mat_expr_genes, :]
                                 .max(axis=1), axis=0)
                         .mean(axis=1))
        res = stats.ttest_rel( mel_in_region, sim_in_region, nan_policy='omit',)
        distro_mels.append(mel_in_region.dropna())
        distro_sims.append(sim_in_region.dropna())

        data['num_{}_change'.format(region_name)] = (
            sum(abs(mel_in_region - sim_in_region)
                > .25))
        data['frac_higher_{}'.format(region_name)] = (
            sum((mel_in_region - sim_in_region) > .25)
            / data['num_{}_change'.format(region_name)]
        )
        data['prob_higher_{}'.format(region_name)] = stats.binom_test(
            [sum((mel_in_region - sim_in_region) > .25),
             sum((sim_in_region - mel_in_region) > .25)]
        )
        print(tf, region, res,
              'mel', mel_in_region.mean(), mel_in_region.ix[tf],
              'sim', sim_in_region.mean(), sim_in_region.ix[tf])

    x_range = np.arange(0, 1, .05)
    if HAS_MPL:
        print("Violinning")
        mpl.figure()
        mpl.violinplot([(dm - ds) for dm, ds in zip(distro_mels, distro_sims)])
        mpl.xticks(np.arange(len(target_regions))+1,
                   [region for region in target_regions],
                   rotation=90)

        mpl.tight_layout()
        mpl.figure()
        ax1 = mpl.gca()
        mpl.figure()
        ax2 = mpl.gca()

        for nmeg, dm, ds, (target, region) in zip(
            nmegs, distro_mels, distro_sims, target_regions.values(),
        ):
            if target == 'hb': continue
            region_diffs = dm - ds
            control_region = (melXmel.ix[target] < 5) & (simXsim.ix[target] < 5)
            control_diffs = ((melXmel.ix[nmeg, control_region]
                              - simXsim.ix[nmeg, control_region])
                             .divide(parental.ix[nmeg, :] .max(axis=1), axis=0)
                            ).mean(axis=1)
            ax1.semilogy(100*x_range,
                         [sum(region_diffs > i)/ sum(region_diffs < -i)
                          for i in x_range],
                         label='{}: {:0.0%} - {:0.0%}'.format(target, *region), basey=2)
            ax2.semilogy(100*x_range,
                         [
                             stats.binom_test([sum(region_diffs > i),
                                               sum(region_diffs < -i)],
                                              p=(sum(control_diffs > i)
                                                 / (sum(control_diffs > i)
                                                    + sum(control_diffs < -i))
                                                )
                                             )
                             for i in x_range
                         ],
                         label='{}: {:0.0%} - {:0.0%}'.format(target, *region),
                         basey=10)
        ax1.legend(loc='upper left')
        ax1.hlines(1, 0, 100)
        ax1.set_ylabel(r'$\#(\Delta > x) \div \#(\Delta < -x)$')
        ax2.hlines(0.05, 0, 100, label='nominal p=.05')
        ax2.legend(loc='lower left')
        ax2.set_ylabel('binomial p value (vs empirical unexpressed control)')

        pu.minimize_ink(ax1)
        pu.minimize_ink(ax2)



    with open('analysis/results/tf_stats.tex', 'w') as outf:
        for var, val in data.items():
            numeric =  isinstance(val, (float, int, int64, int32, int8))
            frac = var.lower().startswith('frac')
            outf.write(create_latex_command(var, val, numeric, frac))



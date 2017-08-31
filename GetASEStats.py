from __future__ import print_function
import numpy as np
import pandas as pd
from collections import Counter
from numpy import isfinite, int64, int32, int16, int8, sign, abs, nan
from scipy.stats import ttest_1samp
import Utils as ut
import PlotUtils as pu
from Utils import startswith, fbgns, pd_kwargs
import HybridUtils as hu
from matplotlib import cm
from progressbar import ProgressBar

def slices_per_embryo(ase):
    return Counter(i.split('_sl')[0] for i in ase.columns)

def create_latex_command(name, value, numeric=False, frac=False):
    name = name.upper().replace('_', '')
    if frac:
        if 0 < abs(value) < .1e-2:
            return '\\newcommand{{\\{}}}[0]{{{:%}}} \n'.format(name, value).replace('%', '\\%')
        return '\\newcommand{{\\{}}}[0]{{{:.1%}}} \n'.format(name, value).replace('%', '\\%').replace('.0', '')
    if numeric:
        return '\\newcommand{{\\{}}}[0]{{{:6,g}}} \n'.format(name, value)
    return '\\newcommand{{\\{}}}[0]{{{}}} \n'.format(name, value)

def get_class(gene, ase, subset='', slices_with_expr=None, expr=None):
    sample = ase.ix[gene]
    sample = sample.select(startswith(subset))

    if slices_with_expr is not None and gene in slices_with_expr.index:
        slices_with_expr = slices_with_expr.ix[gene]
    elif slices_with_expr is None and expr is not None and gene in expr.index:
        slices_with_expr = (expr.ix[gene].select(startswith(subset)) > EXPR_MIN).sum()
    else:
        return nan
    ase_vals = (abs(sample) > ASE_MIN) * sign(sample)
    slices_with_ase = isfinite(sample).sum()
    if slices_with_expr < len(sample) * .90:
        return 99
    if slices_with_ase < .5 * slices_with_expr:
        return 999
    if sum(ase_vals == 1) > slices_with_ase * FRAC_FOR_MATERNAL:
        return 1
    if sum(ase_vals == -1) > slices_with_ase * FRAC_FOR_MATERNAL:
        return -1
    return 0



lott_sort = lambda x: (int(x[1:3]), x[3:])


EXPR_MIN = 10
FRAC_FOR_ASE = 2/3
ASE_MIN = (FRAC_FOR_ASE - (1-FRAC_FOR_ASE))/1
FRAC_FOR_MATERNAL = 0.65

plot_kwargs = {'box_height': 25,
               'col_sep': '_sl',
               'convert': True,
               'draw_box': True,
               'draw_name': True,
               'draw_row_labels': True,
               'make_hyperlinks': True,
               'progress_bar': False,
               'split_columns': True,
               'total_width': 200}

if __name__ == "__main__":
    if 'ase' not in locals() or ('reload_ase' in locals() and locals()['reload_ase']):
        print("Reloading data")
        ase = (pd.read_table('analysis_godot/ase_summary_by_read.tsv', **pd_kwargs)
               .dropna(how='all', axis=1)
               .select(**ut.sel_startswith(('melXsim', 'simXmel')))
              )
        all_ase = ase.copy()
        expr = (pd.read_table('analysis_godot/summary.tsv', **pd_kwargs)
                .drop('---', axis=1, errors='ignore')
                #.dropna(how='all', axis=1)
               )
        lott = pd.read_table('prereqs/journal.pbio.1000590.s002', index_col=0, keep_default_na=False, na_values=[''])
        lott_expr = (lott
                     .ix[:, sorted(lott.columns[5:29], key=lott_sort)]
                     .rename_axis(axis=1, mapper=lambda x: 'lott_sl'+x)
                    )
        reload_ase = False
        to_gn = pd.read_table('prereqs/gene_map_table_fb_2016_01.tsv', index_col=1, skiprows=4).ix[:,0]
        to_fbgn = ut.get_synonyms()

        classes = pd.DataFrame(index=ase.index, columns=['melXsim', 'simXmel'], data=pd.np.nan)
        for ix in ProgressBar()(classes.index):
            for col in classes.columns:
                classes.ix[ix, col] = get_class(ix, ase, subset=col, expr=expr)

    #lott['fbgn'] = to_fbgn[lott.index]
    #lott.drop_duplicates(subset='fbgn', keep=False, inplace=True)
    #lott.index = lott.fbgn

    paris = pd.read_table('prereqs/GSE68062_Gene_CLASS_after_FPKM_normalization.txt', index_col=1,
            header=0,
            names=['gn', 'FBgn', 'yak', 'pse', 'vir']+[s+'_'+n for s in 'mel yak pse vir'.split() for n in 'class prob'.split()],
            )


    in_both = ase.index.intersection(expr.index)
    ase = ase.ix[in_both]
    expr = expr.ix[in_both]
    rn = lambda x: 'parental_' + x.split('_')[2]
    mel_parental = expr.select(**ut.sel_startswith('melXmel')).rename_axis(rn,
                                                                           axis="columns")
    sim_parental = expr.select(**ut.sel_startswith('simXsim')).rename_axis(rn,
                                                                           axis="columns")


    syns = ut.get_synonyms()
    chrom_of = ut.get_chroms(syns)

    males = ('melXsim_cyc14C_rep3', 'simXmel_cyc14C_rep2')
    on_x = [chrom_of[gene] == 'X' for gene in ase.index]
    is_male = [col.startswith(males) for col in ase.columns]
    ase_nomaleX = ase.copy()
    ase_nomaleX.ix[on_x, is_male] = pd.np.nan
    ase = ase_nomaleX


    # Mutliplying your ASE values by parent of origin should make it so that
    # maternal alleles are positive and paternal allels are negative

    parent_of_origin = pd.Series(
            index=ase.columns,
            data=[-1 if c.startswith('m') else 1 for c in ase.columns]
            )


    data = {}
    data['frac_for_ase'] = FRAC_FOR_ASE
    data['frac_for_maternal'] = FRAC_FOR_MATERNAL
    data['expr_min'] = EXPR_MIN

    n_slices = slices_per_embryo(ase)
    data['most_slices'] = max(n_slices.values())
    data['least_slices'] = min(n_slices.values())

    slices_with_expr = (expr > EXPR_MIN).sum(axis=1)
    slices_with_ase = (ase > ASE_MIN).sum(axis=1)
    slices_with_aseval = ase.count(axis=1)
    #slices_with_aseval = slices_with_aseval.where(slices_with_aseval>slices_with_expr, slices_with_expr)
    #slices_with_aseval = slices_with_aseval.where(slices_with_aseval>5, 5)

    bias_dirs = pd.DataFrame(data=pd.np.nan, index=ase.index,
                             columns=['melXsim', 'simXmel'])
    biases = Counter()

    for gene in ase.index:
        gene_biases = hu.is_directionally_biased(ase, gene, two_tailed=True,
                                                 style='cutoff',
                                                 frac_for_biased=FRAC_FOR_MATERNAL,
                                                 ase_level=ASE_MIN,
                                                )
        biases[(gene_biases['melXsim'], gene_biases['simXmel'])] += 1
        bias_dirs.ix[gene] = gene_biases

    maternal = pd.Series(index=ase.index,
                         data=[all(bias_dirs.ix[gene] == (-1, 1)) for gene in
                               ase.index]
                       )
    zygotic = pd.Series(index=ase.index,
                         data=[all(bias_dirs.ix[gene] == (0, 0)) for gene in
                               ase.index]
                       )


    maternal_ttest = pd.Series(index=all_ase.index,
                               data=1.0)
    for gene in maternal_ttest.index:
        if all_ase.ix[gene].count() > 10:
            tstat, pval = ttest_1samp(all_ase.ix[gene] * parent_of_origin, 0,
                                      nan_policy='omit')
            maternal_ttest[gene] = pval# < 1e-1 / len(maternal_ttest)

    print(*ut.true_index(maternal_ttest < (1e-2)), sep='\n',
          file=open('analysis/results/strict_maternal_gene_names.txt', 'w'))
    paternal_mxs = pd.Series(index=ase_nomaleX.index, data=pd.np.nan)
    paternal_sxm = pd.Series(index=ase_nomaleX.index, data=pd.np.nan)
    samedir = pd.Series(index=ase_nomaleX.index, data=pd.np.nan)
    for gene in paternal_sxm.index:
        if expr.ix[gene].max() < EXPR_MIN:
            continue
        gene_ase = ase_nomaleX.ix[gene].multiply(parent_of_origin).dropna()
        per_sample = Counter(col.split('_')[0] for col in gene_ase.index)
        if len(per_sample) < 2 or per_sample.most_common(2)[1][1] < 5:
            continue

        tstat, pval = ttest_1samp(gene_ase.select(startswith('melXsim')), 0)
        paternal_mxs.ix[gene] = pval/2 if tstat < 0 else 1-pval/2
        tstat, pval = ttest_1samp(gene_ase.select(startswith('simXmel')), 0)
        paternal_sxm.ix[gene] = pval/2 if tstat < 0 else 1-pval/2
        #tstat, pval = ttest_ind(gene_ase.select('melXsim'),
                                #gene_ase.select('simXmel'))
    paternal_ttest = paternal_mxs.where(paternal_mxs > paternal_sxm, paternal_sxm).dropna().sort_values(ascending=False)



    data['num_maternal'] = biases[(-1, 1)]
    data['num_paternal'] = biases[(1, -1)]
    data['mel_dominant'] = biases[(-1, -1)]
    data['sim_dominant'] = biases[(1, 1)]

    # Species dominance
    mel_dom = ut.true_index((bias_dirs.T == [-1, -1]).sum() == 2)
    sim_dom = ut.true_index((bias_dirs.T == [ 1,  1]).sum() == 2)

    mel_dom = expr.select(**ut.sel_startswith('simXsim')).T.max()[mel_dom].sort_values().index
    sim_dom = expr.select(**ut.sel_startswith('melXmel')).T.max()[sim_dom].sort_values().index

    print("Making species bias figs")
    pu.svg_heatmap((expr, (sim_parental - mel_parental)/(sim_parental + mel_parental), ase), 'analysis/results/mel_dom.svg',
                   index=mel_dom,
                   norm_rows_by=('maxall', 'center0pre', 'center0pre'),
                   cmap=(pu.ISH, cm.RdBu, cm.RdBu),
                   progress_bar=True,
                   row_labels=[('{:6.1f}'.format(expr.ix[i].max()),
                                chrom_of[i], i)
                               for i in mel_dom],
                   nan_replace='no',
                   **pu.kwargs_heatmap)

    pu.svg_heatmap((expr, (sim_parental - mel_parental)/(sim_parental + mel_parental), ase), 'analysis/results/sim_dom.svg',
                   index=sim_dom,
                   norm_rows_by=('maxall', 'center0pre', 'center0pre'),
                   cmap=(pu.ISH, cm.RdBu, cm.RdBu),
                   progress_bar=True,
                   row_labels=[('{:6.1f}'.format(expr.ix[i].max()),
                                chrom_of[i], i)
                               for i in sim_dom],
                   nan_replace='no',
                   **pu.kwargs_heatmap)


    low_expr_lott = slices_with_expr[lott.index[lott.CLASS == 'mat']] < FRAC_FOR_MATERNAL * len(expr.columns)
    data['lott_maternal_low'] = sum(low_expr_lott)
    has_ase_lott = (slices_with_aseval > FRAC_FOR_MATERNAL * len(ase.columns)).ix[lott.index[lott.CLASS == 'mat']].dropna()
    data['lott_maternal_measured'] = sum(has_ase_lott)

    data['lott_maternal_agree'] = sum(has_ase_lott*maternal.ix[has_ase_lott.index])

    me_mat_lott_zyg = set(lott[lott.CLASS == 'zyg'].index).intersection(maternal[maternal].index)
    me_zyg_lott_mat = set(lott[lott.CLASS == 'mat'].index).intersection(ut.true_index(zygotic))
    data['lott_disagree_t_one'] = len(me_mat_lott_zyg)
    data['lott_disagree_t_two'] = len(me_zyg_lott_mat)
    pu.svg_heatmap((ase.ix[me_mat_lott_zyg], lott_expr.ix[me_mat_lott_zyg]), 'analysis/results/me_mat_lott_zyg.svg',
                   norm_rows_by=('center0pre', 'max'), cmap=(cm.RdBu, cm.viridis),  **plot_kwargs)

    peak_genes = [line.strip() for line in open('analysis/results/peak_genes.txt')]
    logist_genes = [line.strip() for line in open('analysis/results/logist_genes.txt')]

    data['num_peak'] = len(peak_genes)
    data['num_logist'] = len(logist_genes)

    peak_fd = np.fromfile('analysis/results/fd_peak.numpy')
    logist_fd = np.fromfile('analysis/results/fd_logist.numpy')
    peak_r2s = pd.Series.from_csv('analysis/results/all_peak_r2s.csv')
    logist_r2s = pd.Series.from_csv('analysis/results/all_logist_r2s.csv')
    co = 0.45
    data['fd_peak'] = sum(peak_fd > co)
    data['frac_fdr_peak'] = (sum(peak_fd > co) / len(peak_fd)) / (sum(peak_r2s > co) / len(peak_r2s))
    data['frac_max_fdr_peak'] = (1 / len(peak_fd)) / (sum(peak_r2s > co) / len(peak_r2s))
    data['fd_logist'] = sum(logist_fd > co)
    data['frac_fdr_logist'] = sum(logist_fd > co) / len(logist_fd) / (sum(peak_r2s > co) / len(peak_r2s))
    data['frac_max_fdr_logist'] = 1 / len(logist_fd)/ (sum(peak_r2s > co) / len(peak_r2s))

    hb_bind_data = {line.strip()
                    for line in open('analysis/results/hb_wt_emd_0.1.txt')
                    if line.strip() in expr.index}

    #data['hb_diff'] = len(hb_bind_data)
    #data['frac_higher_hb'] = 0

    print(data)

    with open('analysis/results/stats.tex', 'w') as outf:
        for var, val in data.items():
            numeric =  isinstance(val, (float, int, int64, int32, int8))
            frac = var.lower().startswith('frac')
            outf.write(create_latex_command(var, val, numeric, frac))


    if data['num_paternal']:
        pu.svg_heatmap(ase.ix[paternal_ttest.index[paternal_ttest < .05]],
                       'analysis/results/paternal.svg',
                       norm_rows_by='center0pre', cmap=cm.RdBu,
                       hatch_nan=True,hatch_size=1,
                       row_labels=fbgns[paternal_ttest.index],
                       **plot_kwargs)


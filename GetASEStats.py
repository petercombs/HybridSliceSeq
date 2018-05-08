from __future__ import print_function
import numpy as np
import pandas as pd
from collections import Counter
from numpy import isfinite, int64, int32, int16, int8, sign, abs, nan
import Utils as ut
import PlotUtils as pu
from Utils import startswith, fbgns, pd_kwargs
from matplotlib import cm
from warnings import filterwarnings

filterwarnings('ignore', category=FutureWarning)
filterwarnings('ignore', category=DeprecationWarning)

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
               'progress_bar': True,
               'split_columns': True,
               'total_width': 200}

if __name__ == "__main__":
    if 'ase' not in locals() or ('reload_ase' in locals() and locals()['reload_ase']):
        print("Reloading data")
        ase = (pd.read_table('analysis_godot/wasp_summary_by_read.tsv', **pd_kwargs)
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


    in_both = ase.index.intersection(expr.index)
    ase = ase.ix[in_both]
    expr = expr.ix[in_both]
    rn = lambda x: 'parental_' + x.split('_')[2]
    mel_parental = expr.select(**ut.sel_startswith('melXmel')).rename_axis(rn,
                                                                           axis="columns")
    sim_parental = expr.select(**ut.sel_startswith('simXsim')).rename_axis(rn,
                                                                           axis="columns")


    if 'syns' not in locals() or (locals().get('reload_syns', False)):
        syns = ut.get_synonyms()
        chrom_of = ut.get_chroms(syns)
        reload_syns = False

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
    ase_rectified = ase.multiply(parent_of_origin)


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

    print("Species dominance data...")

    deseq_mel = pd.read_table('analysis/mel_deseq.tsv', index_col=0,
                              keep_default_na=False, na_values=['NA', '---'])
    deseq_sim = pd.read_table('analysis/sim_deseq.tsv', index_col=0,
                              keep_default_na=False, na_values=['NA', '---'])

    # Species dominance
    mel_mel_bias = ut.true_index((deseq_mel.padj < .05) & (deseq_mel.log2FoldChange > 0))
    sim_mel_bias = ut.true_index((deseq_sim.padj < .05) & (deseq_sim.log2FoldChange > 0))
    mel_sim_bias = ut.true_index((deseq_mel.padj < .05) & (deseq_mel.log2FoldChange < 0))
    sim_sim_bias = ut.true_index((deseq_sim.padj < .05) & (deseq_sim.log2FoldChange < 0))

    maternal = mel_mel_bias.intersection(sim_sim_bias)
    paternal = mel_sim_bias.intersection(sim_mel_bias)
    mel_dom = mel_mel_bias.intersection(sim_mel_bias)
    sim_dom = mel_sim_bias.intersection(sim_sim_bias)
    zygotic = ut.true_index((deseq_mel.padj > .05) & (deseq_sim.padj > .05))

    data['num_maternal'] = len(maternal)
    data['num_paternal'] = len(paternal)
    data['mel_dominant'] = len(mel_dom)
    data['sim_dominant'] = len(sim_dom)


    mel_dom = expr.select(**ut.sel_startswith('simXsim')).T.max()[mel_dom].sort_values().index
    sim_dom = expr.select(**ut.sel_startswith('melXmel')).T.max()[sim_dom].sort_values().index

    print("Making species bias figs")
    pu.svg_heatmap((expr, (sim_parental - mel_parental)/(sim_parental + mel_parental), ase), 'analysis/results/mel_dom.svg',
                   index=mel_dom,
                   norm_rows_by=('maxall', 'center0pre', 'center0pre'),
                   cmap=(pu.ISH, cm.RdBu, cm.RdBu),
                   progress_bar=True,
                   row_labels=[(
                       '{:6.1f}'.format(
                           expr.ix[i].max() if i in expr.index else np.nan),
                       chrom_of.get(i, '???'),
                       i)
                       for i in mel_dom],
                   nan_replace='no',
                   **pu.kwargs_heatmap)

    pu.svg_heatmap((expr, (sim_parental - mel_parental)/(sim_parental + mel_parental), ase), 'analysis/results/sim_dom.svg',
                   index=sim_dom,
                   norm_rows_by=('maxall', 'center0pre', 'center0pre'),
                   cmap=(pu.ISH, cm.RdBu, cm.RdBu),
                   progress_bar=True,
                   row_labels=[(
                       '{:6.1f}'.format(
                           expr.ix[i].max() if i in expr.index else np.nan),
                       chrom_of.get(i, '???'), i)
                       for i in sim_dom],
                   nan_replace='no',
                   **pu.kwargs_heatmap)


    low_expr_lott = slices_with_expr[lott.index[lott.CLASS == 'mat']] < FRAC_FOR_MATERNAL * len(expr.columns)
    data['lott_maternal_low'] = sum(low_expr_lott)
    has_ase_lott = (slices_with_aseval > FRAC_FOR_MATERNAL * len(ase.columns)).ix[lott.index[lott.CLASS == 'mat']].dropna()
    data['lott_maternal_measured'] = sum(has_ase_lott)

    data['lott_maternal_agree'] = len(maternal.intersection(ut.true_index(lott.CLASS == 'mat')))

    me_mat_lott_zyg = ut.true_index(lott.CLASS == 'zyg').intersection(maternal)
    me_zyg_lott_mat = ut.true_index(lott.CLASS == 'mat').intersection(zygotic)
    me_zyg_lott_mat = ase_rectified.ix[me_zyg_lott_mat].T.mean().sort_values().index
    data['lott_disagree_t_one'] = len(me_mat_lott_zyg)
    data['lott_disagree_t_two'] = len(me_zyg_lott_mat)
    #pu.svg_heatmap((ase.ix[me_mat_lott_zyg], lott_expr.ix[me_mat_lott_zyg]),
    #               'analysis/results/me_mat_lott_zyg.svg',
    #               norm_rows_by=('center0pre', 'max'),
    #               cmap=(cm.RdBu, cm.viridis),
    #               **plot_kwargs)

    small_heatmap_kwargs = plot_kwargs.copy()
    small_heatmap_kwargs['box_height'] = 5
    small_heatmap_kwargs['draw_row_labels'] = False
    pu.svg_heatmap((ase.ix[me_mat_lott_zyg], lott_expr.ix[me_mat_lott_zyg]),
                   'analysis/results/me_mat_lott_zyg.svg',
                   norm_rows_by=('center0pre', 'max'),
                   cmap=(cm.RdBu, cm.viridis),
                   **plot_kwargs)
    pu.svg_heatmap((ase.ix[me_zyg_lott_mat], lott_expr.ix[me_zyg_lott_mat]),
                   'analysis/results/me_zyg_lott_mat.svg',
                   norm_rows_by=('center0pre', 'max'),
                   cmap=(cm.RdBu, cm.viridis),
                   **small_heatmap_kwargs)

    peak_genes = [line.strip() for line in open('analysis/results/asepeak_genes.txt')]
    logist_genes = [line.strip() for line in open('analysis/results/aselogist_genes.txt')]

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

    print(data)

    with open('analysis/results/stats.tex', 'w') as outf:
        for var, val in data.items():
            numeric =  isinstance(val, (float, int, int64, int32, int16, int8))
            frac = var.lower().startswith('frac')
            outf.write(create_latex_command(var, val, numeric, frac))


    if data['num_paternal']:
        pu.svg_heatmap(ase.ix[paternal],
                       'analysis/results/paternal.svg',
                       norm_rows_by='center0pre', cmap=cm.RdBu,
                       hatch_nan=True,hatch_size=1,
                       row_labels=fbgns[paternal],
                       **plot_kwargs)


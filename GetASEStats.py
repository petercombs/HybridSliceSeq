import pandas as pd
from collections import Counter
from numpy import isfinite, int64, int32, int16, int8, sign, abs, nan
from scipy.stats import ttest_1samp
import Utils
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
        return '\\newcommand{{\\{}}}[0]{{{:.1%}}} \n'.format(name, value).replace('%', '\\%').replace('.0', '')
    if numeric:
        return '\\newcommand{{\\{}}}[0]{{{:,}}} \n'.format(name, value)
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
               .select(**Utils.sel_startswith(('melXsim', 'simXmel')))
              )
        expr = pd.read_table('analysis_godot/summary_fb.tsv', **pd_kwargs).dropna(how='all', axis=1)
        lott = pd.read_table('prereqs/journal.pbio.1000590.s002', index_col=0, keep_default_na=False, na_values=[''])
        reload_ase = False
        to_gn = pd.read_table('prereqs/gene_map_table_fb_2016_01.tsv', index_col=1, skiprows=4).ix[:,0]
        to_fbgn = Utils.get_synonyms()

        lott['fbgn'] = to_fbgn[lott.index]
        lott.drop_duplicates(subset='fbgn', keep=False, inplace=True)
        lott.index = lott.fbgn

        paris = pd.read_table('prereqs/GSE68062_Gene_CLASS_after_FPKM_normalization.txt', index_col=1,
                header=0,
                names=['gn', 'FBgn', 'yak', 'pse', 'vir']+[s+'_'+n for s in 'mel yak pse vir'.split() for n in 'class prob'.split()],
                )

        classes = pd.DataFrame(index=ase.index, columns=['melXsim', 'simXmel'], data=pd.np.nan)
        for ix in ProgressBar()(classes.index):
            for col in classes.columns:
                classes.ix[ix, col] = get_class(ix, ase, subset=col, expr=expr)

    in_both = ase.index.intersection(expr.index)
    ase = ase.ix[in_both]
    expr = expr.ix[in_both]

    chrom_of = {}
    for row in open('prereqs/gene_map_table_fb_2016_01.tsv'):
        if row.startswith('#') or not row.strip():
            continue
        data = row.split()
        chrom_of[data[1]] = data[-1].split(':')[0]

    males = ('melXsim_cyc14C_rep3', 'simXmel_cyc14C_rep2')
    on_x = [chrom_of[gene] == 'X' for gene in ase.index]
    is_male = [col.startswith(males) for col in ase.columns]
    ase_nomaleX = ase.copy()
    ase_nomaleX.ix[on_x, is_male] = pd.np.nan


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
    slices_with_aseval = isfinite(ase).sum(axis=1)
    slices_with_aseval = slices_with_aseval.where(slices_with_aseval>slices_with_expr, slices_with_expr)
    slices_with_aseval = slices_with_aseval.where(slices_with_aseval>5, 5)

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

    low_expr_lott = slices_with_expr[lott.index[lott.CLASS == 'mat']] < FRAC_FOR_MATERNAL * len(expr.columns)
    data['lott_maternal_low'] = sum(low_expr_lott)
    has_ase_lott = (slices_with_aseval > FRAC_FOR_MATERNAL * len(ase.columns)).ix[lott.index[lott.CLASS == 'mat']].dropna()
    data['lott_maternal_measured'] = sum(has_ase_lott)

    data['lott_maternal_agree'] = sum(has_ase_lott*maternal.ix[has_ase_lott.index])

    mat_lott_zyg = set(lott[lott.CLASS == 'zyg'].index).intersection(maternal[maternal].index)
    data['lott_disagree'] = len(mat_lott_zyg)



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

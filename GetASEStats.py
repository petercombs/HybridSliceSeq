import pandas as pd
from collections import Counter, defaultdict
from numpy import isfinite, int64, int32, int16, int8, sign, abs, nan
import Utils
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
    sample = sample.select(Utils.startswith(subset))

    if slices_with_expr is not None and gene in slices_with_expr.index:
        slices_with_expr = slices_with_expr.ix[gene]
    elif slices_with_expr is None and expr is not None and gene in expr.index:
        slices_with_expr = (expr.ix[gene].select(Utils.startswith(subset)) > EXPR_MIN).sum()
    else:
        return nan
    ase_vals = (abs(sample) > ASE_MIN) * sign(sample)
    if slices_with_expr < len(sample) * .90:
        return 99
    if sum(ase_vals == 1) > slices_with_expr * FRAC_FOR_MATERNAL:
        return 1
    if sum(ase_vals == -1) > slices_with_expr * FRAC_FOR_MATERNAL:
        return -1
    return 0



    


EXPR_MIN = 10
FRAC_FOR_ASE = 2/3
ASE_MIN = (FRAC_FOR_ASE - (1-FRAC_FOR_ASE))/1
FRAC_FOR_MATERNAL = 0.65

if __name__ == "__main__":
    kwargs = dict(index_col=0,
            keep_default_na=False, 
            na_values=['---'])


    if 'ase' not in locals() or ('reload_ase' in locals() and locals()['reload_ase']):
        print("Reloading data")
        ase = pd.read_table('analysis_godot/ase_summary.tsv', **kwargs).dropna(how='all', axis=1)
        expr = pd.read_table('analysis_godot/summary_fb.tsv', **kwargs).dropna(how='all', axis=1)
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

    maternal = (ase.multiply(parent_of_origin) > ASE_MIN).sum(axis=1) > (FRAC_FOR_MATERNAL * slices_with_aseval)
    paternal = (ase.multiply(parent_of_origin) < -ASE_MIN).sum(axis=1) > (FRAC_FOR_MATERNAL * slices_with_aseval)

    data['num_maternal'] = sum(maternal)
    data['num_paternal'] = sum(paternal)

    low_expr_lott = slices_with_expr[lott.index[lott.CLASS == 'mat']] < FRAC_FOR_MATERNAL * len(expr.columns)
    data['lott_maternal_low'] = sum(low_expr_lott)
    has_ase_lott = (slices_with_aseval > FRAC_FOR_MATERNAL * len(ase.columns)).ix[lott.index[lott.CLASS == 'mat']].dropna()
    data['lott_maternal_measured'] = sum(has_ase_lott)

    data['lott_maternal_agree'] = sum(has_ase_lott*maternal.ix[has_ase_lott.index])

    mat_lott_zyg = set(lott[lott.CLASS == 'zyg'].index).intersection(maternal[maternal].index)
    data['lott_disagree'] = len(mat_lott_zyg)


    mel_dominant = (ase < -ASE_MIN).sum(axis=1) > FRAC_FOR_MATERNAL * slices_with_aseval
    data['mel_dominant'] = sum(mel_dominant)
    sim_dominant = (ase > ASE_MIN).sum(axis=1) > FRAC_FOR_MATERNAL * slices_with_aseval
    data['sim_dominant'] = sum(sim_dominant)

    print(data)

    with open('analysis/results/stats.tex', 'w') as outf:
        for var, val in data.items():
            numeric =  isinstance(val, (float, int, int64, int32, int8))
            frac = var.lower().startswith('frac')
            outf.write(create_latex_command(var, val, numeric, frac))

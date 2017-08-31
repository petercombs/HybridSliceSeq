import pandas as pd
from Utils import sel_startswith, pd_kwargs, get_xs, get_chroms
from FitASEFuncs import (fit_all_ase, logistic, peak,
                         calculate_variance_explained)
from multiprocessing import Pool
import PlotUtils as pu
from matplotlib import cm
from progressbar import ProgressBar as pb

male_hybrid_embryos = ('melXsim_cyc14C_rep3', 'simXmel_cyc14C_rep2')
female_hybrid_embryos = ('melXsim_cyc14C_rep1', 'melXsim_cyc14C_rep2',
                         'simXmel_cyc14C_rep1')

if __name__ == "__main__":
    expr = pd.read_table('godot/summary_fb.tsv', **pd_kwargs)
    ase = (pd.read_table('godot/ase_summary.tsv', **pd_kwargs)
           .dropna(how='all', axis=0)
          )

    chrom_of = get_chroms()

    ase = ase.select(lambda x: chrom_of[x] != 'X')

    expr_males = expr.select(**sel_startswith(male_hybrid_embryos))
    expr_females = expr.select(**sel_startswith(female_hybrid_embryos))

    ase_males = ase.select(**sel_startswith(male_hybrid_embryos))
    ase_females = ase.select(**sel_startswith(female_hybrid_embryos))

    ase_xs = get_xs(ase)
    ase_maternals = pd.Series(
        index=ase_xs.index,
        data=[1 if col.startswith('simXmel') else -1 for col in ase_xs.index]
    )

    if 'logistic_females' in locals() and locals().get('recalculate', True):
        with Pool() as p:
            logistic_females = fit_all_ase(ase_females, logistic,
                                           ase_xs.ix[ase_females.columns],
                                           pool=p, progress=True)
            peak_females = fit_all_ase(ase_females, peak,
                                       ase_xs.ix[ase_females.columns],
                                       pool=p, progress=True)

        female_logistic_r2 = calculate_variance_explained(
            ase_females, ase_xs.ix[ase_females.columns],
            logistic,
            logistic_females,
        )
        female_peak_r2 = calculate_variance_explained(
            ase_females, ase_xs.ix[ase_females.columns],
            peak,
            peak_females
        )

        male_logistic_r2 = calculate_variance_explained(
            ase_males, ase_xs.ix[ase_males.columns],
            logistic,
            logistic_females,
        ).clip(0, 1)
        male_peak_r2 = calculate_variance_explained(
            ase_males, ase_xs.ix[ase_males.columns],
            peak,
            peak_females
        ).clip(0, 1)
        recalculate = False

    pu_kwargs = {
        'box_height': 60,
        'col_sep': '_sl',
        'convert': True,
        'draw_box': True,
        'draw_name': False,
        'draw_row_labels': True,
        'make_hyperlinks': True,
        'max_width': 880,
        'progress_bar': False,
        'split_columns': True,
        'total_width': 200,
        'nan_replace' : 0.5,
        'vspacer': 0}

    diffs = set()
    for gene, diff in (female_logistic_r2 - male_logistic_r2).items():
        if female_logistic_r2[gene] > .4 and abs(diff) > .25:
            diffs.add(gene)
    for gene, diff in (female_peak_r2 - male_peak_r2).items():
        if female_peak_r2[gene] > .4 and abs(diff) > .25:
            diffs.add(gene)

    for gene in pb()(diffs):
        pu.svg_heatmap(
            (
                None, expr_females.ix[[gene]],
                None, ase_females.ix[gene],
                None, ase_males.select(**sel_startswith('melXsim')).ix[gene],
                None, ase_males.select(**sel_startswith('simXmel')).ix[gene],
                None, expr_males.select(**sel_startswith('melXsim')).ix[[gene]],
                None, expr_males.select(**sel_startswith('simXmel')).ix[[gene]],
            ),
            'analysis_godot/results/sex_diff/{}.svg'.format(gene),
            norm_rows_by=(
                'female expression', 'max',
                'females - L{:.03f} P{:.03f}'.format(female_logistic_r2[gene],
                                                     female_peak_r2[gene]),
                'center0pre',
                'males - L{:.03f} P{:.03f}'.format(male_logistic_r2[gene],
                                                   male_peak_r2[gene]),
                'center0pre', '', 'center0pre',
                'male expression', 'max', '', 'max',
            ),
            cmap=(
                None, pu.ISH,
                None, cm.RdBu,
                None, cm.RdBu, None, cm.RdBu,
                None, pu.ISH, None, pu.ISH,
            ),
            **pu_kwargs
        )




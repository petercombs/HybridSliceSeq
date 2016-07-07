import pandas as pd
from Utils import sel_startswith, pd_kwargs, get_xs
from FitASEFuncs import (fit_all_ase, logistic, peak,
                         calculate_variance_explained)
from multiprocessing import Pool

male_hybrid_embryos = ('melXsim_cyc14C_rep3', 'simXmel_cyc14C_rep2')
female_hybrid_embryos = ('melXsim_cyc14C_rep1', 'melXsim_cyc14C_rep2',
                         'simXmel_cyc14C_rep1')

if __name__ == "__main__":
    expr = pd.read_table('godot/summary_fb.tsv', **pd_kwargs)
    ase = (pd.read_table('godot/ase_summary.tsv', **pd_kwargs)
           .dropna(how='all', axis=0)
          )

    ase_males = ase.select(**sel_startswith(male_hybrid_embryos))
    ase_females = ase.select(**sel_startswith(female_hybrid_embryos))

    ase_xs = get_xs(ase)
    ase_maternals = pd.Series(
        index=ase_xs.index,
        data=[1 if col.startswith('simXmel') else -1 for col in ase_xs.index]
    )

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



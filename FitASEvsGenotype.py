from __future__ import print_function
from numpy import (exp, array,  nan, inf)
from multiprocessing import Pool
from Utils import (sel_startswith, get_xs, pd_kwargs, get_synonyms,
                   startswith, get_chroms)
import pandas as pd
import numpy as np
import PlotUtils as pu
from FitASEFuncs import (logistic, peak, fit_all_ase,
                         calculate_variance_explained)

if __name__ == "__main__":
    expr = pd.read_table('analysis_godot/summary.tsv', **pd_kwargs).dropna(how='all', axis=1)
    ase = (pd
           .read_table('analysis_godot/ase_summary_by_read.tsv',
                       **pd_kwargs
                       )
           .dropna(how='all', axis=1)
           .dropna(how='all', axis=0)
           .select(**sel_startswith(('melXsim', 'simXmel')))
          )
    ase_limited = ase.select(**sel_startswith('melXsim'))
    chrom_of = get_chroms()

    males = ('melXsim_cyc14C_rep3', 'simXmel_cyc14C_rep2')
    on_x = chrom_of[ase.index] == 'X'
    is_male = [col.startswith(males) for col in ase.columns]
    ase.ix[on_x, is_male] = np.nan


    xs = get_xs(ase)
    xs_ltd = get_xs(ase_limited)
    colnames = ['Amp', 'width', 'center', 'y_offset']
    recalc_ase = locals().get('recalc_ase', True)
    if recalc_ase:
        with Pool() as p:
            res_logist = fit_all_ase(ase, logistic, xs, colnames, p,
                                     progress=True).dropna()
            res_logist_limited = fit_all_ase(ase_limited, logistic, xs_ltd, colnames, p,
                                             progress=True).dropna()


            res_peak = fit_all_ase(ase, peak, xs, colnames, p,
                                   progress=True).dropna()
            res_peak_limited = fit_all_ase(ase_limited, peak, xs_ltd, colnames, p,
                                           progress=True).dropna()


    recalc_ase = False



    r2_logist = calculate_variance_explained(ase, xs, logistic, res_logist)
    r2_logist_perm = calculate_variance_explained(ase_limited, xs, logistic,
                                                  res_logist_limited)
    r2_peak = calculate_variance_explained(ase, xs, peak, res_peak)
    r2_peak_perm = calculate_variance_explained(ase_limited, xs, peak,
                                                res_peak_limited)



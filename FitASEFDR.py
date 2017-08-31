import numpy as np
import pandas as pd
from FitASEFuncs import (logistic, peak, fit_all_ase,
                         calculate_variance_explained)

from multiprocessing import Pool
from progressbar import ProgressBar as pbar
from Utils import (sel_startswith, get_xs, pd_kwargs, get_chroms)



if __name__ == "__main__":
    expr = pd.read_table('analysis_godot/summary_fb.tsv', **pd_kwargs).dropna(how='all', axis=1)
    ase = (pd
           .read_table('analysis_godot/ase_summary_by_read.tsv',
                       **pd_kwargs
                       )
           .dropna(how='all', axis=1)
           .dropna(how='all', axis=0)
           .select(**sel_startswith(('melXsim', 'simXmel')))
          )
    chrom_of = get_chroms()

    males = ('melXsim_cyc14C_rep3', 'simXmel_cyc14C_rep2')
    on_x = [chrom_of[gene] == 'X' if gene in chrom_of else False for gene in ase.index]
    is_male = [col.startswith(males) for col in ase.columns]
    ase.ix[on_x, is_male] = np.nan


    xs = get_xs(ase)
    colnames = ['Amp', 'width', 'center', 'y_offset']
    peak_r2s = []
    logist_r2s = []

    with Pool() as p:
        for i in pbar(max_value=1000)(range(1000)):
            new_xs = pd.Series(index=xs.index,
                               data=np.random.permutation(xs))
            res_logist = fit_all_ase(ase, logistic, new_xs, colnames, p,
                                     progress=False).dropna()
            res_peak = fit_all_ase(ase, peak, new_xs, colnames, p,
                                   progress=False).dropna()

            peak_r2s.extend(calculate_variance_explained(
                ase, new_xs, peak, res_peak
            ))

            logist_r2s.extend(calculate_variance_explained(
                ase, new_xs, logistic, res_logist
            ))

    np.array(peak_r2s).tofile('analysis/results/fd_peak.numpy')
    np.array(logist_r2s).tofile('analysis/results/fd_logist.numpy')

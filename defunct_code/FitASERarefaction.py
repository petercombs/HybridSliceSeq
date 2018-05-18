from numpy import (isfinite)
from scipy.stats import linregress
from multiprocessing import Pool
from collections import defaultdict
from glob import glob
from os import path
import pandas as pd
import numpy as np
from progressbar import ProgressBar
from FitASE import (logistic, peak, fit_all_ase,
                    calculate_variance_explained)


if __name__ == "__main__":
    ase = (pd
           .read_table('analysis_godot/ase_summary.tsv',
                       index_col=0,
                       keep_default_na=False, na_values=['---'],)
           .dropna(how='all', axis=1)
          )
    max_slice = defaultdict(int)
    for sl in ase.columns:
        sl = sl.split('_')
        emb = '_'.join(sl[:2])
        max_slice[emb] = max(max_slice[emb], int(sl[2][2:]))

    colnames = ['Amp', 'width', 'center', 'y_offset']
    all_peaks = {}
    all_logis = {}
    with Pool() as p:
        for fname in glob('analysis_godot/ase_summary_by_read_in_*subset.tsv'):
            print(fname)
            ase = (pd
                   .read_table(fname,
                               index_col=0,
                               keep_default_na=False,
                               na_values=['---', '-', ''])
                   .dropna(how='all', axis=1)
                   .dropna(how='all', axis=0)
                  )
            print(ase.shape)
            xs = pd.Series(index=ase.columns,
                           data=[int(a.split('_')[2][2:])/max_slice['_'.join(a.split('_')[:2])] for a in ase.columns if 'sl' in a])
            res_logist = fit_all_ase(ase, logistic, xs, colnames, p).dropna()
            res_peak = fit_all_ase(ase, peak, xs, colnames, p).dropna()

            res_lin = pd.DataFrame(
                index=ase.index,
                columns=['slope', 'intercept', 'pval', 'r'],
                data=np.nan
                )

            for gene in ProgressBar()(ase.index):
                cols = isfinite(ase.ix[gene])
                linreg = linregress(xs[cols], ase.ix[gene, cols])
                if linreg.pvalue < .05:
                    res_lin.ix[gene] = [linreg.slope, linreg.intercept, linreg.pvalue, linreg.rvalue]



            r2_logist = calculate_variance_explained(ase, xs, logistic, res_logist)
            r2_peak = calculate_variance_explained(ase, xs, peak, res_peak)
            all_peaks[path.basename(fname)] = r2_peak.index[r2_peak>0.5]
            all_logis[path.basename(fname)] = r2_logist.index[r2_logist>.5]


    print("Found {: 4} {:<10} of which {} are likely false"
          .format(sum(r2_logist > 0.5), 'logistic', None))
    print("Found {: 4} {:<10} of which {} are likely false"
          .format(sum(r2_peak > 0.5), 'peak', None))

    good_amps_logist = res_logist.Amp[r2_logist>0.5].sort_values(inplace=False)
    good_amps_peak = res_peak.Amp[r2_peak>0.5].sort_values(inplace=False)
    #good_amps_logist.to_csv('analysis/results/logist.tsv', sep='\t')
    #good_amps_peak.to_csv('analysis/results/peak.tsv', sep='\t')
    res_logist.ix[r2_logist > 0.5].to_csv('analysis/results/logist.tsv', sep='\t')
    res_peak.ix[r2_peak > 0.5].to_csv('analysis/results/peak.tsv', sep='\t')






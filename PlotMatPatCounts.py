from __future__ import print_function
from glob import glob
from collections import Counter
from os import path
import pandas as pd
import Utils as ut
from scipy import stats
from matplotlib.pyplot import (scatter, xlabel, ylabel, gca, savefig, plot,
                               xticks, yticks, tight_layout)
import warnings
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=DeprecationWarning)

if __name__ == "__main__":
    dirs = glob('analysis_godot/*sl*')
    dirshort = [path.basename(d) for d in dirs]
    counts = pd.DataFrame(index=sorted(dirshort), columns=[None, 0, -1, 1])
    hybs = ('melXsim', 'simXmel')

    curr = Counter()
    for dirname in dirs:
        l = next(open(path.join(dirname, 'wasp_gene_ase_by_read.tsv')))
        exec('curr = ' + l.split(maxsplit=4)[4].strip())
        counts.ix[path.basename(dirname)] = curr

    mat_counts = [counts.ix[d, (1 if d.startswith('sim') else -1)] for d in counts.index]
    pat_counts = [counts.ix[d, (1 if d.startswith('mel') else -1)] for d in counts.index]
    mat_counts = pd.Series(index=counts.index, data=mat_counts)
    pat_counts = pd.Series(index=counts.index, data=pat_counts)

    colors = dict(melXsim='r', simXmel='b', melXmel='gray', simXsim='gray')
    scatter(mat_counts, pat_counts, c=[colors[d.split('_')[0]] for d in counts.index])
    counts.to_csv('analysis/results/allele_counts.tsv', sep='\t')
    print('Paternal Counts', pat_counts.select(ut.startswith(hybs)).mean())

    xlabel('Maternal Counts')
    ylabel('Paternal Counts')
    ax = gca()
    fit = stats.linregress(mat_counts.select(ut.startswith(hybs)), pat_counts.select(ut.startswith(hybs)))
    xmin, xmax = mat_counts.min(), mat_counts.select(ut.startswith(hybs)).max()
    plot([xmin, xmax], [xmin * fit.slope + fit.intercept, xmax * fit.slope + fit.intercept], 'r:')
    ax = gca()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_bounds(pat_counts.min(), pat_counts.max())
    ax.spines['bottom'].set_bounds(mat_counts.min(), mat_counts.max())
    xt = [mat_counts.min(), 200000, 400000, 600000, 800000, mat_counts.max()]
    yt = [pat_counts.select(ut.startswith(hybs)).min(), 40000, 80000, 120000,
          pat_counts.max()]
    xticks(xt, ['{:,}'.format(int(i)) for i in xt])
    yticks(yt, ['{:,}'.format(int(i)) for i in yt])
    tight_layout()

    savefig('analysis/results/matpat_counts.png', dpi=300)

from __future__ import print_function
from numpy import (exp, array, isfinite, nan, inf)
from scipy.optimize import curve_fit
from scipy.stats import linregress, ttest_1samp
from multiprocessing import Pool
from collections import Counter
from Utils import (sel_startswith, get_xs, pd_kwargs, fbgns, get_synonyms,
                   startswith, get_chroms)
import pandas as pd
import itertools as it
import numpy as np
import inspect
from progressbar import (ProgressBar, Percentage, Timer, Bar, FileTransferSpeed,
                        SimpleProgress, AdaptiveETA)
import PlotUtils as pu
import Utils as ut
import warnings
if __name__ == "__main__":
    from matplotlib import cm
    import matplotlib.pyplot as mpl

def logistic(x, A, k, x0, y0):
    x = array(x)
    return A/(1+exp(k*(x-x0))) - y0

def peak(x, A, w, x0, y0):
    x = array(x)
    return A * exp(-(x-x0)**2/w**2) - y0

def deriv_peak(x, A, w, x0, y0):
    return - A * (x-x0) * exp(-(x-x0)**2/w**2) / w**2

def fit_func(func, index, data, xs, p0=None, median_in=None, randomize=False,
             print_error=False):
    warnings.filterwarnings('ignore', '.*overflow encountered.*')
    warnings.filterwarnings('ignore', '.*Covariance of the parameters.*')
    ys = data.ix[index]
    if median_in and not median_in[0] < ys.median() < median_in[1]:
        print("Median {} outside of range [{}, {}]"
              .format(ys.median(), *median_in)
             )
        return array([nan, nan, nan, nan])
    if randomize:
        xs = np.random.permutation(xs)
    keep = array(isfinite(ys))
    if sum(keep) < len(xs)/3:
        if print_error:
            print("Not enough values to fit ASE ({}/{})"
                 .format(sum(keep), len(xs)))
        return array([nan, nan, nan, nan])
    if p0 is None:
        pos_amp = ys.max() - ys.median()
        neg_amp = ys.min() - ys.median()
        if pos_amp > -neg_amp:
            amp = pos_amp
        else:
            amp = neg_amp
        p0 = [amp, 20, xs.mean(), ys.median()]
    w_min = 1
    w_max = 200
    x0_min = -0.1
    x0_max = 1.1
    if func == peak:
        w_min = 0.15
        w_max = 0.6
        x0_min = 0.0
        x0_max = 1.0
        p0[1] = 0.4
        if (ys.max() - ys.median()) < (ys.median() - ys.min()):
            p0[0] = ys.median() - ys.min()
            p0[2] = xs[ys.idxmin()]
        else:
            p0[2] = xs[ys.idxmax()]
    else:
        pass
    try:
        if print_error:
            print(p0)
        res = curve_fit(func,
                         array(xs[keep]),
                         array(ys[keep]),
                         p0,
                         bounds=[[-2, w_min, x0_min, -1],
                                 [2, w_max, x0_max, 1]],
                        )[0]
        for low, high, val in zip([-2, w_min, x0_min, -1],
                                  [2, w_max, x0_max, 1],
                                  res):
            if not low <= val <= high:
                return array([nan, nan, nan, nan])
            return res

    except (ValueError, RuntimeError) as err:
        if print_error:
            print(err)
        return array([nan, nan, nan, nan])

def fit_all_ase(ase, func, xs, colnames=None, pool=None, progress=False,
                median_in=None):
    if colnames is None:
        colnames = inspect.getargspec(func).args
        colnames.remove('x')
    if pool is not None:
        if isinstance(pool, int):
            pool = Pool(processes=pool)
        if progress:
            results = [pool.apply_async(fit_func, (func, i, ase, xs),
                                        {'median_in': median_in},)
                       for i in ase.index]
            res = []
            pbar = ProgressBar(maxval=len(results), widgets=[ "|",
                                                             Percentage(), "|",
                                                             SimpleProgress(),
                                                             Bar(),
                                                             FileTransferSpeed(), "|",
                                                             Timer(), "|",
                                                             AdaptiveETA(),
                                                            ])
            for i, r in zip(ase.index, pbar(results)):
                res.append(r.get())
        else:
            res = list(pool.starmap(
                fit_func,
                zip(
                    it.repeat(func),
                    ase.index,
                    it.repeat(ase),
                    it.repeat(xs),
                    it.repeat(None), # p0
                    it.repeat(median_in),
                ),
            ))

    else:
        res = list(it.starmap(
            fit_func,
            zip(
                it.repeat(func),
                ase.index,
                it.repeat(ase),
                it.repeat(xs),
                it.repeat(None), # p0
                it.repeat(median_in),
                )
            ))

    return pd.DataFrame(
            index=ase.index,
            columns=colnames,
            data=res
            )


def calculate_variance_explained(ase, xs, func, params):
    r2 = pd.Series(index=params.index, data=np.inf)
    for ix in params.index:
        r2.ix[ix] = 1-(
            ((ase.ix[ix] - func(xs, *params.ix[ix]))**2).sum()
            / ((ase.ix[ix] - ase.ix[ix].mean())**2).sum()
        )

    return r2


overlap_genes = ['CG43355',  # with sala
                 'CG45087',  # with Pepck
                 'CG30080',  # with CG42662
                ]
# Genes with significantly overlapping exons, but different CDSes. Called
# manually based on identical ASE for the genes.

males = ('melXsim_cyc14C_rep3', 'simXmel_cyc14C_rep2')

def parse_args():
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('--expression-file', '-x',
                        default='analysis_godot/summary.tsv')
    parser.add_argument('--ase-file', '-a',
                        default='analysis_godot/ase_summary_by_read.tsv')
    parser.add_argument('--prefix', '-p', default='')
    parser.add_argument('--suffix', '-s', default='')
    parser.add_argument('--overlapping-genes', default=overlap_genes, nargs='*')
    parser.add_argument('--male-samples', default=males, nargs='+')
    parser.add_argument('--print-keggs', default=False, action='store_true')
    parser.add_argument('--cutoff-r2', default=0.45, type=float)
    parser.add_argument('--min-samples', default=20, type=int)
    parser.add_argument('--min-var', default=False, type=float)
    parser.add_argument('--fit-ymin', default=-1, type=float)
    parser.add_argument('--fit-ymax', default=1, type=float)
    args = parser.parse_args()
    args.male_samples = tuple(args.male_samples)
    return args



if __name__ == "__main__":
    args = parse_args()
    expr = pd.read_table(args.expression_file, **pd_kwargs).drop('---', axis=1, errors='ignore')
    ase = (pd
           .read_table(args.ase_file,
                       **pd_kwargs
                       )
           .dropna(how='all', axis=1)
           .dropna(how='all', axis=0)
           .drop(args.overlapping_genes, errors='ignore')
           .select(**sel_startswith(('melXsim', 'simXmel')))
          )
    ase_perm = pd.DataFrame(
            data=np.random.permutation(ase.T).T,
            index=ase.index,
            columns=ase.columns,
            )
    chrom_of = get_chroms()

    if args.male_samples and 'keep' not in args.male_samples:
        on_x = chrom_of[ase.index] == 'X'
        is_male = [col.startswith(args.male_samples) for col in ase.columns]
        ase.ix[on_x, is_male] = np.nan
    ase = ase.ix[ase.T.count() >= args.min_samples]
    if args.min_var:
        ase = ase.ix[ase.T.var() >= args.min_var]
    ase_perm = ase_perm.ix[ase.index]


    xs = get_xs(ase)
    colnames = ['Amp', 'width', 'center', 'y_offset']
    recalc_ase = locals().get('recalc_ase', True)
    if recalc_ase:
        with Pool() as p:
            res_logist = fit_all_ase(ase, logistic, xs, colnames, p,
                                     progress=True).dropna()
            res_logist_perm = pd.DataFrame(0, index=res_logist.index,
                                           columns=res_logist.columns)
            #res_logist_perm = fit_all_ase(ase_perm, logistic, xs, colnames, p,
                                          #progress=True).dropna()
            res_peak = fit_all_ase(ase, peak, xs, colnames, p,
                                   progress=True).dropna()
            res_peak_perm = pd.DataFrame(0, index=res_peak.index,
                                           columns=res_peak.columns)
            #res_peak_perm = fit_all_ase(ase_perm, peak, xs, colnames, p,
                                        #progress=True).dropna()

            res_lin = pd.DataFrame(
                    index=ase.index,
                    columns=['slope', 'intercept', 'pval', 'r'],
                    data=np.nan
                    )
            res_lin_perm = res_lin.copy()

            pbar = ProgressBar()
            for gene in pbar(ase.index):
                cols = isfinite(ase.ix[gene])
                if sum(cols) == 0:
                    continue
                linreg = linregress(xs[cols], ase.ix[gene, cols])
                if linreg.pvalue < .05:
                    res_lin.ix[gene] = [linreg.slope, linreg.intercept, linreg.pvalue, linreg.rvalue]

                #cols = isfinite(ase_perm.ix[gene])
                #linreg = linregress(xs[cols], ase_perm.ix[gene, cols])
                #if linreg.pvalue < .05:
                    #res_lin_perm.ix[gene] = [linreg.slope, linreg.intercept, linreg.pvalue, linreg.rvalue]
    recalc_ase = False



    r2_logist = calculate_variance_explained(ase, xs, logistic, res_logist)
    r2_logist_perm = np.zeros_like(r2_logist)
    r2_peak = calculate_variance_explained(ase, xs, peak, res_peak)
    r2_peak_perm = np.zeros_like(r2_peak)


    co = args.cutoff_r2

    print("Found {: 4} {:<10} of which {} are likely false".format(sum(r2_logist > co), 'logistic', sum(r2_logist_perm > co)))
    print("Found {: 4} {:<10} of which {} are likely false".format(sum(r2_peak > co), 'peak', sum(r2_peak_perm > co)))

    good_amps_logist = res_logist.Amp[r2_logist>co].sort_values(inplace=False)
    good_amps_peak = res_peak.Amp[r2_peak>co].sort_values(inplace=False)
    #good_amps_logist.to_csv('analysis/results/logist.tsv', sep='\t')
    #good_amps_peak.to_csv('analysis/results/peak.tsv', sep='\t')
    res_logist.ix[r2_logist > co].to_csv(
        'analysis/results/{prefix}logist{suffix}.tsv'
        .format(prefix=args.prefix, suffix=args.suffix),
        sep='\t')
    res_peak.ix[r2_peak > co].to_csv(
        'analysis/results/{prefix}peak{suffix}.tsv'
        .format(prefix=args.prefix, suffix=args.suffix),
        sep='\t')

    kwargs = dict(
        progress_bar=False,
        col_sep='_sl',
        total_width=150,
        box_height=25,
        split_columns=True,
        draw_box=True,
        draw_row_labels=True,
        draw_name=True,
        make_hyperlinks=True,
        convert=True,
    )

    r2_peak_genes = {
        gene:r2_peak.ix[gene] for gene in r2_peak.index
        if ((r2_peak.ix[gene] > co)
            and (gene not in r2_logist.index
                 or not (r2_peak.ix[[gene]] < r2_logist.ix[[gene]]).all()))
    }
    r2_peak_genes = res_peak.Amp[r2_peak_genes].sort_values().index

    r2_logist_genes = {
        gene:r2_logist.ix[gene] for gene in r2_logist.index
        if ((r2_logist.ix[gene] > co)
            and not (r2_logist.ix[[gene]] < r2_peak.ix[[gene]]).all())
    }
    r2_logist_genes = res_logist.Amp[r2_logist_genes].sort_values().index

    print('\n'.join(r2_peak_genes),
          file=open('analysis/results/{prefix}peak{suffix}_genes.txt'
                    .format(prefix=args.prefix, suffix=args.suffix),
                    'w'))
    print('\n'.join(r2_logist_genes),
          file=open('analysis/results/{prefix}logist{suffix}_genes.txt'
                    .format(prefix=args.prefix, suffix=args.suffix),
                    'w'))

    best_r2 = r2_peak.copy()
    for gene in r2_logist.index:
        # Note that I use "and not > " rather than "<" to catch nans
        if ((gene in best_r2 and not best_r2[gene] > r2_logist[gene])
            or (gene not in best_r2)):
            best_r2[gene] = r2_logist[gene]

    best_r2.to_csv('analysis/results/{prefix}svase{suffix}_best'.format(
        prefix=args.prefix, suffix=args.suffix), sep='\t')


    pu.svg_heatmap(ase.ix[r2_logist_genes],
                   'analysis/results/{prefix}logist{suffix}_ase.svg'
                   .format(prefix=args.prefix, suffix=args.suffix),
                   norm_rows_by='center0pre',
                   row_labels=['{:.02f} {}'.format(best_r2[g], g) for g in
                               r2_logist_genes],
                   cmap=cm.RdBu,
                   **kwargs)

    pu.svg_heatmap(ase.ix[r2_peak_genes],
                   'analysis/results/{prefix}peak{suffix}_ase.svg'
                   .format(prefix=args.prefix, suffix=args.suffix),
                   norm_rows_by='center0pre',
                   row_labels=['{:.02f} {}'.format(best_r2[g], g) for g in
                               r2_peak_genes],
                   cmap=cm.RdBu,
                   **kwargs)

    pu.svg_heatmap(ase.ix[r2_logist_genes].select(**ut.sel_contains('rep1')),
                   'analysis/results/{prefix}logist{suffix}_ase_r1.svg'
        .format(prefix=args.prefix, suffix=args.suffix),
                   norm_rows_by='center0pre',
                   cmap=cm.RdBu,
                   **kwargs)

    pu.svg_heatmap(ase.ix[r2_peak_genes].select(**ut.sel_contains('rep1')),
                   'analysis/results/{prefix}peak{suffix}_ase_r1.svg'
                   .format(prefix=args.prefix, suffix=args.suffix),
                   norm_rows_by='center0pre',
                   cmap=cm.RdBu,
                   **kwargs)
    pu.svg_heatmap(expr.ix[r2_logist_genes].select(**sel_startswith(('melXsim',
                                                                    'simXmel'))),
                   'analysis/results/{prefix}logist{suffix}_expr_hyb.svg'
                   .format(prefix=args.prefix, suffix=args.suffix),
                   norm_rows_by='maxall',
                   **kwargs)

    pu.svg_heatmap(expr.ix[r2_peak_genes].select(**sel_startswith(('melXsim',
                                                                  'simXmel'))),
                   'analysis/results/{prefix}peak{suffix}_expr_hyb.svg'
                   .format(prefix=args.prefix, suffix=args.suffix),
                   norm_rows_by='maxall',
                   **kwargs)

    pu.svg_heatmap(expr.ix[r2_logist_genes],
                   'analysis/results/{prefix}logist{suffix}_expr.svg'
                   .format(prefix=args.prefix, suffix=args.suffix),
                   norm_rows_by='maxall',
                   **kwargs)

    pu.svg_heatmap(expr.ix[r2_peak_genes],
                   'analysis/results/{prefix}peak{suffix}_expr.svg'
                   .format(prefix=args.prefix, suffix=args.suffix),
                   norm_rows_by='maxall',
                   **kwargs)

    mpl.figure(figsize=(2, len(r2_logist_genes)))
    for i, gene in enumerate(r2_logist_genes):
        mpl.subplot(len(r2_logist_genes), 1, i+1)
        mpl.plot(sorted(xs), logistic(array(sorted(xs)), *res_logist.ix[gene]),
             'g-', linewidth=3)
        mpl.yticks([])
        mpl.xticks([])
        mpl.xlim(0, 1)
        mpl.ylim(args.fit_ymin, args.fit_ymax)
    mpl.tight_layout(pad=0, h_pad=0, w_pad=0)
    mpl.savefig('analysis/results/{prefix}logist{suffix}_fits'
                .format(prefix=args.prefix, suffix=args.suffix),)

    mpl.figure(figsize=(2, len(r2_peak_genes)))
    for i, gene in enumerate(r2_peak_genes):
        mpl.subplot(len(r2_peak_genes), 1, i+1)
        mpl.plot(sorted(xs), peak(array(sorted(xs)), *res_peak.ix[gene]),
             'g-', linewidth=3)
        mpl.yticks([])
        mpl.xticks([])
        mpl.xlim(0, 1)
        mpl.ylim(args.fit_ymin, args.fit_ymax)
    mpl.tight_layout(pad=0, h_pad=0, w_pad=0)
    mpl.savefig('analysis/results/{prefix}peak{suffix}_fits'
                .format(prefix=args.prefix, suffix=args.suffix),
               )

    if args.print_keggs:
        synonyms = get_synonyms()
        wnt_genes = [line.strip() for line in open('prereqs/wnt.kegg.genes')]
        wnt_scores = pd.Series(index=synonyms[wnt_genes],
                               data=best_r2[synonyms[wnt_genes]])
        wnt_scores.index = ['dme:Dmel_' + CG for CG in wnt_genes]
        wnt_scores.index.name = '#dme'
        wnt_scores.name = 'svASE'
        (wnt_scores
         .sort_values(na_position='first')
         .to_csv('analysis/results/wnt_scores.tsv',
                 sep='\t',
                 header=True)
        )

        all_cgs = synonyms.select(startswith(('CG1', 'CG2', 'CG3', 'CG4', 'CG5',
                                              'CG6', 'CG7', 'CG8', 'CG9')))
        all_scores = pd.Series(index=all_cgs,
                               data=best_r2[all_cgs])
        all_scores.index = ['dme:Dmel_' + CG for CG in all_cgs.index]
        all_scores.index.name = '#dme'
        all_scores.name = 'svASE'
        all_scores2 = all_scores.copy()
        all_scores2.index = [ix.split(':')[1] for ix in all_scores2.index]
        (all_scores
         .sort_values(na_position='first')
         .dropna()
         .to_csv('analysis/results/all_svase_scores_cg.tsv',
                 sep='\t',
                 header=True)
        )
        keggs = {line.split()[0]:line.split()[1].strip().strip(',').split(',')
                 for line in open('prereqs/kegg_database.txt')}
        kegg_stats = Counter()
        kegg_pvals = Counter()
        for pathway in ProgressBar()(keggs):
            ttest_res = ttest_1samp(all_scores2[keggs[pathway]].dropna(),
                                    0.045130029332897552)
            if len(all_scores2[keggs[pathway]].dropna()) > 5:
                kegg_stats[pathway] = ttest_res.statistic
                kegg_pvals[pathway] = ttest_res.pvalue
                if ttest_res.statistic > 0:
                    pathway_gene_names = synonyms[[ix.split('_')[1]
                                                   for ix in keggs[pathway]]]
                    pathway_svases = best_r2[pathway_gene_names].sort_values(ascending=False)
                    pathway_name = pathway.split('~')[1].replace('_/', '')
                    pu.svg_heatmap(ase.ix[pathway_svases.index],
                                   'analysis/results/{}_ase.svg'.format(pathway_name),
                                   row_labels=fbgns.ix[pathway_svases.index],
                                   norm_rows_by='center0pre',
                                   cmap=cm.RdBu,
                                   **kwargs
                                  )

                    pu.svg_heatmap(expr.ix[pathway_svases.index],
                                   'analysis/results/{}_expr.svg'.format(pathway_name),
                                   row_labels=fbgns.ix[pathway_svases.index],
                                   norm_rows_by='max',
                                   cmap=pu.ISH,
                                   **kwargs
                                  )

"""A script to find trans-regulatory differences causing svASE

The basic idea here is that if, between the two different directions of the
cross, there are radically different svASE functions, there have to be
trans-regulatory changes that are causing it.

So, what I'll be doing is fitting functions to each cross separately, then see
if there's a difference by... not sure yet. Maybe see if the R^2 is lower? How
much lower? Maybe I'll just sort by the difference in the variance explained and
eyeball it.

"""


import pandas as pd
from multiprocessing import Pool
from collections import defaultdict
from Utils import sel_startswith
from FitASEFuncs import (calculate_variance_explained, logistic, peak,
                         fit_all_ase)
import PlotUtils as pu
import progressbar as pb
from matplotlib import cm

kwargs = dict(
    norm_rows_by='center0pre',
    progress_bar=False,
    col_sep='_sl',
    total_width=200,
    box_height=60,
    split_columns=True,
    draw_box=True,
    draw_row_labels=True,
    draw_name=False,
    make_hyperlinks=True,
    convert=True,
    vspacer=0,
    max_width=200,
    cmap=cm.RdBu,
)


if __name__ == "__main__":
    ase = (pd
           .read_table('analysis_godot/ase_summary_by_read.tsv',
                       index_col=0,
                       keep_default_na=False, na_values=['---'],)
           .dropna(how='all', axis=1)
           .dropna(how='all', axis=0)
           .select(**sel_startswith(('melXsim', 'simXmel')))
          )
    paris = pd.read_table('prereqs/GSE68062_Gene_CLASS_after_FPKM_normalization.txt',
                  index_col=1)['mel.CLASS']
    pzyg = paris[paris == 'zyg']

    melXsim = ase.select(**sel_startswith('melXsim')).select(pzyg.__contains__)
    simXmel = ase.select(**sel_startswith('simXmel')).select(pzyg.__contains__)

    fbgns = pd.read_table('prereqs/gene_map_table_fb_2016_01.tsv',
                          index_col=1,skiprows=5).ix[:, 0]

    max_slice = defaultdict(int)
    for sl in ase.columns:
        sl = sl.split('_sl')
        emb = sl[0]
        max_slice[emb] = max(max_slice[emb], int(sl[1][0:2]))

    xs = pd.Series(index=ase.columns,
                   data=[int(a.split('_sl')[1][:2])/max_slice[a.split('_sl')[0]]
                         for a in ase.columns if 'sl' in a])
    xs_mxs = xs[melXsim.columns]
    xs_sxm = xs[simXmel.columns]

    colnames = ['Amp', 'width', 'center', 'y_offset']
    with Pool() as p:
        res_logist_mxs = fit_all_ase(melXsim, logistic, xs_mxs,
                                     colnames, p,
                                     progress=True,
                                     median_in=(-.8, .8),
                                    ).dropna()
        print("Fit 1/4 done")

        r2_logist_mxs_native = calculate_variance_explained(
            melXsim, xs_mxs, logistic, res_logist_mxs)
        r2_logist_mxs_other  = calculate_variance_explained(
            simXmel, xs_sxm, logistic, res_logist_mxs)



        res_peak_mxs = fit_all_ase(melXsim, peak, xs_mxs,
                                   colnames, p,
                                   progress=True,
                                   median_in=(-.8, .8),
                                  ).dropna()
        print("Fit 2/4 done")

        r2_peak_mxs_native = calculate_variance_explained(melXsim, xs_mxs,
                                                     peak, res_peak_mxs)
        r2_peak_mxs_other  = calculate_variance_explained(simXmel, xs_sxm,
                                                     peak, res_peak_mxs)

        res_logist_sxm = fit_all_ase(simXmel, logistic, xs_sxm,
                                     colnames, p,
                                     progress=True,
                                     median_in=(-.8, .8),
                                    ).dropna()
        print("Fit 3/4 done")

        r2_logist_sxm_native = calculate_variance_explained(
            simXmel, xs_sxm, logistic, res_logist_sxm)
        r2_logist_sxm_other  = calculate_variance_explained(
            melXsim, xs_mxs, logistic, res_logist_sxm)



        res_peak_sxm = fit_all_ase(simXmel, peak, xs_sxm,
                                   colnames, p,
                                   progress=True,
                                   median_in=(-.8, .8),
                                  ).dropna()
        print("Fit 4/4 done")

        r2_peak_sxm_native = calculate_variance_explained(
            simXmel, xs_sxm, peak, res_peak_sxm)
        r2_peak_sxm_other  = calculate_variance_explained(
            melXsim, xs_mxs, peak, res_peak_sxm)

    r2_peak_diff = (r2_peak_mxs_native - r2_peak_mxs_other.clip(0,1)).dropna().sort_values()
    r2_logist_diff = (r2_logist_mxs_native -
                      r2_logist_mxs_other.clip(0,1)).dropna().sort_values()

    pbar = pb.ProgressBar()
    for gene in pbar(r2_peak_diff.select(pzyg.__contains__).index[-50:]):
        val = r2_peak_diff.ix[gene]
        pu.svg_heatmap(ase.ix[gene],
                       'analysis/results/transdiff/pd_mxs_{:.02f}-{}.svg'.format(val,
                                                                         fbgns[gene]),
                       **kwargs
                      )

    pbar.finish()
    pbar = pb.ProgressBar()
    for gene in pbar(r2_logist_diff.select(pzyg.__contains__).index[-50:]):
        val = r2_logist_diff.ix[gene]
        pu.svg_heatmap(ase.ix[gene],
                       'analysis/results/transdiff/ld_mxs_{:.02f}-{}.svg'.format(val,
                                                                         fbgns[gene]),
                       **kwargs
                      )


    r2_sxm_peak_diff = (r2_peak_sxm_native - r2_peak_sxm_other.clip(0,1)).dropna().sort_values()
    r2_sxm_logist_diff = (r2_logist_sxm_native -
                      r2_logist_sxm_other.clip(0,1)).dropna().sort_values()

    pbar = pb.ProgressBar()
    for gene in pbar(r2_sxm_peak_diff.select(pzyg.__contains__).index[-50:]):
        val = r2_sxm_peak_diff.ix[gene]
        pu.svg_heatmap(ase.ix[gene],
                       'analysis/results/transdiff/pd_sxm_{:.02f}-{}.svg'.format(val,
                                                                         fbgns[gene]),
                       **kwargs
                      )

    pbar.finish()
    pbar = pb.ProgressBar()
    for gene in pbar(r2_sxm_logist_diff.select(pzyg.__contains__).index[-50:]):
        val = r2_sxm_logist_diff.ix[gene]
        pu.svg_heatmap(ase.ix[gene],
                       'analysis/results/transdiff/ld_sxm_{:.02f}-{}.svg'.format(val,
                                                                         fbgns[gene]),
                       **kwargs
                      )

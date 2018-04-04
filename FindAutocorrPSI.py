from tqdm import tqdm
import pandas as pd
from numpy import arange, nan, mean
from collections import defaultdict
import numpy as np
import Utils as ut
import PlotUtils as pu
import CluToGene as spliceid
import multiprocessing as mp
import matplotlib.cm as cm
from os import getcwd
import re

number = re.compile("[0-9]+")

def parse_fasta_data(fasta_fname):
    transcripts_by_gene = defaultdict(list)
    transcript_lens = {}
    for line in open(fasta_fname):
        if not line.startswith('>'): continue
        fbtr, *annot = line.strip().split(' ')
        fbtr = fbtr.strip('>')
        annots = dict(a.strip().strip(';').split('=')
                      for a in annot
                      if '=' in a)
        first, *rest, last = number.findall(annots['loc'].split(':')[1])
        transcript_lens[fbtr] = abs((int(last) - int(first) ))+1
        transcripts_by_gene[annots['parent']].append(fbtr)
    return pd.Series(transcript_lens), transcripts_by_gene

def estimate_pvals(psi, ase, n_reps, min_periods=None):
    #try:
        if min_periods is None:
            min_periods = len(psi)/5
        in_both = np.isfinite(psi*ase)
        if sum(in_both) < min_periods:
            return nan
        psi = psi[in_both]
        ase = ase[in_both]
        observed = psi.corr(ase, min_periods=min_periods)
        if not np.isfinite(observed):
            return nan
        random = []
        for i in range(n_reps):
            np.random.shuffle(psi)
            random.append(psi.corr(ase, min_periods=min_periods))
        random = pd.Series(list(sorted(np.abs(random))))
        return (n_reps-random.searchsorted(abs(observed), 'right')[0])/n_reps
    #except Exception as e:
        #print(psi.name)
        #raise e

def estimate_all_pvals(psi, ase, n_reps, min_periods=None, pool=None,
                       progress=True):
    if pool is None:
        pool = mp.Pool()
    if progress:
        pbar = tqdm
    else:
        pbar = lambda x: x
    if 'estimate_pvals' not in locals():
        from FindAutocorrPSI import estimate_pvals

    jobs = {}
    assert np.all(psi.columns == ase.columns)
    res = pd.Series(index=psi.index, data=np.nan)
    for ix in pbar(psi.index):
        jobs[ix] = pool.apply_async(estimate_pvals,
                                    (psi.ix[ix],
                                     ase.ix[[ut.fbgns[ix.split('_')[0].split('+')[0]]]].squeeze(),
                                     n_reps, min_periods))
    for ix in pbar(psi.index):
        res[ix] = jobs[ix].get()
    return res

cluster_args = dict(time= '2:30:00', mem='40G',
                    partition='owners,hns,normal,hbfraser',
                    scriptpath='logs', outpath='logs', runpath=getcwd(),
                    cpus=4, cores=4)

def fyrd_estimate_pvals(psi, ase, n_reps, min_periods=None,
                        n_genes_per_job=100):
    import fyrd
    outs = {}
    jobs = []
    for i in range(0, len(psi), n_genes_per_job):
        jobs.append(fyrd.submit(estimate_all_pvals,
                                (psi.iloc[i:i+n_genes_per_job],
                                 ase, n_reps, min_periods),
                               **cluster_args))

    for i in tqdm(list(range(len(jobs)))):
        job = jobs[i]
        res = job.get()
        for ix in res.index:
            outs[ix] = res[ix]

    return pd.Series(outs)

def transcripts_from_exon(exons_gtf):
    retval = {}
    for line in open(exons_gtf):
        data = line.strip().split('\t')
        if data[2] != 'exonic_part': continue
        annots = ut.parse_annotation(data[-1])
        exon_id = '{}_{}'.format(annots['gene_id'],
                                 annots['exonic_part_number'])
        retval[exon_id] = annots['transcripts'].split('+')
    return retval


def dist_from_exon_to_transcript_end(reference_gtf, exons_gtf, progress=False):
    if progress:
        pbar = tqdm
    else:
        pbar = lambda x: x
    total_transcript_size = defaultdict(lambda : (1e10, -1))
    transcripts_per_gene = defaultdict(list)

    for line in pbar(open(reference_gtf)):
        if line.startswith('#'): continue
        data = line.strip().split('\t')
        if data[2] != 'exon': continue
        annots = ut.parse_annotation(data[-1])
        FBtr = annots['transcript_id']
        FBgn = annots['gene_id']
        left, right = total_transcript_size[FBtr]
        left = min(left, int(data[3]))
        right = max(left, int(data[4]))
        total_transcript_size[FBtr] = (left, right)
        transcripts_per_gene[FBgn].append(FBtr)

    for transcript in pbar(total_transcript_size):
        left, right = total_transcript_size[transcript]
        total_transcript_size[transcript] = (right - left)
        optional_transcript_length = {}

    for gene in pbar(transcripts_per_gene):
            transcripts = transcripts_per_gene[gene]
            tlens = [total_transcript_size[t]
                     for t in transcripts]
            minlen = min(tlens)
            maxlen = max(tlens)
            for transcript in transcripts:
                optional_transcript_length[transcript] = (
                    (total_transcript_size[transcript] - minlen)
                    / (maxlen - minlen + 1))

    exons = defaultdict(list)

    for line in pbar(open(exons_gtf)):
        data = line.strip().split('\t')
        if data[2] != 'exonic_part': continue
        annots = ut.parse_annotation(data[-1])
        exon_id = "{}_{}".format(annots['gene_id'], annots['exonic_part_number'])
        transcripts = annots['transcripts'].split('+')
        for transcript in transcripts:
            exons[exon_id].append(optional_transcript_length[transcript])
    return exons

if __name__ == "__main__":
    if 'chrom_of' not in locals():
        chrom_of = ut.get_chroms()

    txlens, txs_by_gene = parse_fasta_data('prereqs/dmel-all-transcript-r5.57.fasta')

    psi = (pd.read_table('analysis_godot/psi_summary.tsv',
                              **ut.pd_kwargs)
           .select(**ut.sel_contains(('melXsim', 'simXmel'))))
    ase = (pd.read_table('analysis_godot/ase_summary_by_read.tsv',
                        **ut.pd_kwargs)
           .select(**ut.sel_startswith(('melXsim', 'simXmel'))))
    on_x = chrom_of[ase.index] == 'X'
    is_male = [col.startswith(('melXsim_cyc14C_rep3', 'simXmel_cyc14C_rep2')) for col in ase.columns]
    ase.ix[on_x, is_male] = np.nan

    rectified_ase = ase.multiply([1 if ix.startswith('melXsim') else -1
                                  for ix in ase.columns])


    rectified_ase.columns = psi.columns

    pv10k=fyrd_estimate_pvals(psi, rectified_ase, 10000,
                                              n_genes_per_job=500)
    pv10k.to_csv('analysis/results/pv10k.csv')

#    rectified_ase_by_exons = pd.DataFrame(index=psi.index, columns=psi.columns,
#                                         data=np.nan)
#    for fb in tqdm(ut.fbgns.index):
#        gn = ut.fbgns[fb]
#        if gn not in rectified_ase.index: continue
#        starts_with = rectified_ase_by_exons.index.map(ut.startswith(fb))
#        rectified_ase_by_exons.loc[starts_with, :] = rectified_ase.ix[gn]

    if 'ac_many' not in locals():
        max_ac = 10
        xs = ut.get_xs(psi)
        sxs = xs.sort_values()
        psi_x_sorted = psi.loc[:, sxs.index]
        ac_many = pd.DataFrame(index=psi.index, columns=arange(1,max_ac),
                               data={i:
                                     [psi_x_sorted.loc[gene].dropna().autocorr(i)
                                      for gene in tqdm(psi.index)] for i in
                                     arange(1,max_ac)})

        psi_rand = psi.copy()
        for ix in tqdm(psi_rand.index):
            is_good = np.isfinite(psi_rand.ix[ix])
            dat = np.array(psi_rand.loc[ix, is_good])
            np.random.shuffle(dat)
            psi_rand.ix[ix, is_good] = dat
        ac_many_rand = pd.DataFrame(index=psi.index, columns=arange(1,max_ac),
                                    data={i:
                                          [psi_rand.loc[gene].dropna().autocorr(i)
                                           for gene in tqdm(psi.index)] for
                                          i in arange(1,max_ac)})


    psi_counts = psi.T.count() > psi.shape[1]/3
    plist = ut.true_index(ac_many.loc[psi_counts].T.mean() > .146)
    zyg_corrs = pd.Series(index=plist,
                          data=[psi.loc[ex].corr(
                              rectified_ase.loc[spliceid.get_genes_in_exon(ex).split('_')[0].split('+')[0]],
                              min_periods=30
                          )
                              for ex in plist ])
                          #data=[psi.loc[ex].corr(rectified_ase.ix[[ex]]) for ex
                                #in plist])
    plist = zyg_corrs.sort_values().index

    geneset = {g for gs in plist for g in gs.split('_')[0].split('+')}

    if 'optional_exon_lens' not in locals():
        optional_exon_lens = dist_from_exon_to_transcript_end('Reference/mel_good.gtf',
                                         'Reference/mel_good_exons.gtf', True)
        optional_exon_lens = pd.Series(optional_exon_lens)

    pu.svg_heatmap((
        ase.ix[[spliceid.get_genes_in_exon(ex).split('_')[0].split('+')[0]
                            for ex in plist]],
        #ase.ix[plist],
        psi.ix[plist]),
                   'analysis/results/psi-autocorr-fdr5.svg',
                   cmap=(cm.RdBu, cm.viridis),
                   norm_rows_by=('center0pre', 'fullvar'),
                   row_labels=[('{:.03f}'.format(ac_many.loc[ex].mean()),
                                psi.loc[ex].min(),
                                psi.loc[ex].max(),
                                #ex.split('_')[1],
                                '{:.1f}kb'.format(min(txlens[txs_by_gene[ex.split('_')[0].split('+')[0]]])/1000),
                                '{:.1f}kb'.format(max(txlens[txs_by_gene[ex.split('_')[0].split('+')[0]]])/1000),
                                #'{:.02f}'.format(mean(optional_exon_lens[ex])),
                                '{:.02f}'.format(zyg_corrs[ex]),
                                #'{:.02f}'.format(psi.loc[ex].corr(rectified_ase.loc[spliceid.get_genes_in_exon(ex).split('_')[0].split('+')[0]])),
                                spliceid.get_genes_in_exon(ex))
                               for ex in
                               plist],
                   **pu.kwargs_heatmap)

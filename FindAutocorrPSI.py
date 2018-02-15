from tqdm import tqdm
import pandas as pd
from numpy import arange, nan
from collections import defaultdict
import numpy as np
import Utils as ut
import PlotUtils as pu
import CluToGene as spliceid

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

    rectified_ase_by_exons = pd.DataFrame(index=psi.index, columns=psi.columns,
                                         data=np.nan)
    for fb in tqdm(ut.fbgns.index):
        gn = ut.fbgns[fb]
        if gn not in rectified_ase.index: continue
        starts_with = rectified_ase_by_exons.index.map(ut.startswith(fb))
        rectified_ase_by_exons.loc[starts_with, :] = rectified_ase.ix[gn]

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
            np.shuffle(dat)
            psi_rand.ix[ix, is_good] = dat
        ac_many_rand = pd.DataFrame(index=psi.index, columns=arange(1,max_ac),
                                    data={i:
                                          [psi_rand.loc[gene].dropna().autocorr(i)
                                           for gene in tqdm.tqdm(psi.index)] for
                                          i in arange(1,max_ac)})


    psi_counts = psi.T.count() > psi.shape[1]/3
    plist = ut.true_index(ac_many.loc[psi_counts].T.mean() > .146)
    zyg_corrs = pd.Series(index=plist,
                          #data=[psi.loc[ex].corr(rectified_ase.loc[spliceid.get_genes_in_exon(ex).split('_')[0].split('+')[0]])
                                #for ex in plist ])
                          data=[psi.loc[ex].corr(rectified_ase.ix[[ex]]) for ex
                                in plist])
    plist = zyg_corrs.sort_values().index

    geneset = {g for gs in plist for g in gs.split('_')[0].split('+')}

    pu.svg_heatmap((
        #ase.ix[[spliceid.get_genes_in_exon(ex).split('_')[0].split('+')[0]
                            #for ex in plist]],
        ase.ix[plist],
        psi.ix[plist]),
                   'analysis/results/psi-autocorr-fdr5.svg',
                   cmap=(cm.RdBu, cm.viridis),
                   norm_rows_by=('center0pre', 'fullvar'),
                   row_labels=[('{:.03f}'.format(ac_many.loc[ex, 1]),
                                psi.loc[ex].min(),
                                psi.loc[ex].max(),
                                ex.split('_')[1],
                                '{:.02f}'.format(min(optional_exon_lens[ex])),
                                '{:.02f}'.format(max(optional_exon_lens[ex])),
                                '{:.02f}'.format(mean(optional_exon_lens[ex])),
                                '{:.02f}'.format(zyg_corrs[ex]),
                                #'{:.02f}'.format(psi.loc[ex].corr(rectified_ase.loc[spliceid.get_genes_in_exon(ex).split('_')[0].split('+')[0]])),
                                spliceid.get_genes_in_exon(ex))
                               for ex in
                               plist],
                   **pu.kwargs_heatmap)

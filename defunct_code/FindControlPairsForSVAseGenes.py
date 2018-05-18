from __future__ import print_function
import Utils as ut
import pandas as pd
from scipy import stats
import DistributionDifference as dd
from progressbar import ProgressBar as pbar
from multiprocessing import Pool

def get_ortholog_TSS_data(target='Dsim'):
    mel_tss = {}
    target_tss = {}
    strands = {'1': 0, '-1': 1}
    for line in open('prereqs/gene_orthologs_fb_2015_03.tsv'):
        if not line.strip() or line.startswith('#'): continue
        data = line.split()
        gene, chrom, pos, strand = data[1:5]
        target_gene, target_chrom, target_pos, target_strand = data[6:10]
        if target_gene.startswith(target):
            pos = int(pos.split('..')[strands[strand]])
            target_pos = int(target_pos.split('..')[strands[target_strand]])
            mel_tss[gene] = ('dmel_'+chrom, pos)
            target_tss[gene] = (target.lower()+'_'+target_chrom, target_pos)
    return mel_tss, target_tss

if __name__ == "__main__":
    ase = (pd.read_table('analysis_godot/ase_summary_by_read.tsv', **ut.pd_kwargs)
           .select(**ut.sel_startswith(('melXsim', 'simXmel'))))
    best_r2s = pd.Series.from_csv('analysis/results/svase_best', sep='\t')
    expr = pd.read_table('godot/summary.tsv', **ut.pd_kwargs).drop('---', axis=1)
    mel_tss, sim_tss = get_ortholog_TSS_data()

    has_svase = (ut.true_index(best_r2s.sort_values(ascending=False) > .25)
                 .intersection(mel_tss.keys()))
    no_svase = (best_r2s.index[best_r2s < .01]
                .intersection(expr.index)
                .intersection(mel_tss.keys()))
    median_expr = expr.T.median()

    num_sim_expr = pd.Series(index=has_svase, data=-1)
    already_used = set()
    best_match = pd.DataFrame(index=has_svase, data={'gene': '', 'emd': 1.0})
    p = Pool()
    for gene in pbar()(has_svase):
        similar_expr = ut.true_index((.5 * median_expr[gene] < median_expr[no_svase])
                                     & (median_expr[no_svase] < 2 * median_expr[gene]))
        similar_expr = similar_expr.difference(already_used)
        diff_jobs = {
            target: p.apply_async(dd.earth_mover_multi,
                                  (expr.loc[gene], expr.loc[target]))
            for target in similar_expr
        }
        pattern_diffs = pd.Series({target: diff_jobs[target].get()
                                   for target in similar_expr}).sort_values()
        best_match.loc[gene, 'gene'] = pattern_diffs.index[0]
        best_match.loc[gene, 'emd'] = pattern_diffs[0]
        already_used.add(pattern_diffs.index[0])

    best_match.index.name='svase_gene'
    best_match.to_csv('analysis/results/transposon_dists/best_matches.tsv',
                      sep='\t')
    svase_tss_file = open('analysis/results/transposon_dists/svase_tss.bed', 'w')
    nosvase_tss_file = open('analysis/results/transposon_dists/nosvase_tss.bed', 'w')
    for gene in ut.true_index(best_match.emd < .1):
        print(mel_tss[gene][0], mel_tss[gene][1], mel_tss[gene][1]+1, gene,
              sep='\t', file=svase_tss_file)
        print(sim_tss[gene][0].replace('Scf_', ''),
              sim_tss[gene][1], sim_tss[gene][1]+1, gene,
              sep='\t', file=svase_tss_file)
        ogene = best_match.gene[gene]
        print(mel_tss[ogene][0], mel_tss[ogene][1],
              mel_tss[ogene][1]+1, ogene+'<--'+gene,
              sep= '\t', file=nosvase_tss_file)
        print(sim_tss[ogene][0].replace('Scf_', ''),
              sim_tss[ogene][1], sim_tss[ogene][1]+1,
              ogene+'<--'+gene,
              sep= '\t', file=nosvase_tss_file)

    svase_tss_file.close()
    nosvase_tss_file.close()


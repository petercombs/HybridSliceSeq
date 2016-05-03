from pysam import Samfile
import pandas as pd
from sys import argv
from os import path
import numpy as np
import editdistance as ed
from collections import defaultdict

isfinite = lambda x: x is not np.nan

complement = {
        ord('A'): ord('T'),
        ord('G'): ord('C'),
        ord('C'): ord('G'),
        ord('T'): ord('A'),
        }

sample = 'MA50'

def get_matching_snps(bamfile, chrom, sorted_snps):
    snps_matching = pd.DataFrame(columns=sorted_snps)
    for read in bamfile.fetch(chrom, snps.iloc[0], snps.iloc[-1]):
        aln_seq = read.query_alignment_sequence
        if read.is_read2:
            aln_seq = aln_seq.translate(complement)
        for pos, base in zip(read.get_reference_positions(), aln_seq):
            if pos+1 in snps_matching.columns:
                snps_matching.ix[read.qname, pos+1] = base
    return snps_matching

def get_linked_snps(snps_matching, n_linked=2):
    num_snps = snps_matching.applymap(isfinite).sum(axis=1)
    return snps_matching.ix[num_snps >= n_linked].dropna(how='all', axis=1)

def class_refalt(snps_in_read, ref_alt_dict, seen_snps = set()):
    snps_in_read = snps_in_read.dropna()
    snp_pos = snps_in_read.index
    snps_in_read = ''.join(snps_in_read)
    altseq = ''.join(ref_alt_dict[i-1][1] for i in snp_pos)
    if (snps_in_read == altseq) or (snps_in_read == altseq.translate(complement)):
        return 1
    refseq = ''.join(ref_alt_dict[i-1][0] for i in snp_pos)
    if (snps_in_read == refseq) or (snps_in_read == refseq.translate(complement)):
        return -1
    if len(snps_in_read) < 3:
        pass
    elif seen_snps is not None and snps_in_read not in seen_snps:
        diff_ed = ed.eval(snps_in_read, refseq) - ed.eval(snps_in_read, altseq)
        print(snps_in_read, refseq, altseq, ['==', '=>', '<='][np.sign(diff_ed)], sep='\t')
        seen_snps.add(snps_in_read)
    return 0
    


if __name__ == "__main__":
    snps = pd.read_table('../on_mel/true_hets.tsv', header=None, names=['chrom', 'pos'])
    ref_alt = defaultdict(dict)
    for var_line in open('../on_mel/melsim_variant.bed'):
        chrom, pos, _, refalt = var_line.split()
        ref_alt[chrom][int(pos)] = refalt.split('|')
        


    chr_input = argv[1]
    pos_low = argv[2]
    pos_hi = argv[3]

    snps = snps.query('chrom == "{}" and {} <= pos <= {}'.format(chr_input, pos_low, pos_hi)).pos.sort_values()
    melXsim = get_matching_snps(Samfile('../MA50/assigned_dmelR.bam'), chr_input, snps)
    simXmel = get_matching_snps(Samfile('../RA24/assigned_dmelR.bam'), chr_input, snps)
    mel = get_matching_snps(Samfile('../on_mel/mel_gdna_bowtie2.bam'), chr_input, snps)
    sim = get_matching_snps(Samfile('../on_mel/sim_gdna_bowtie2.bam'), chr_input, snps)







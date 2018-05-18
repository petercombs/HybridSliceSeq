import pandas as pd
from progressbar import ProgressBar as pb
from collections import defaultdict
from sys import stdout

dist = 1e4

def get_snps_and_len(tss, dnase, snps, seen_dnase, exons):
    num_snps = 0
    len_regions = 0
    crms = dnase.query('chrom == "{chr}" and ((abs(start - {pos}) < {dist}) or (abs(stop-{pos}) < {dist}))'
            .format(chr=tss.chr.replace('dmel_', ''),
                pos=tss.TSS_start,
                dist=dist))
    for i, crm in crms.iterrows():
        if i in seen_dnase: continue
        seen_dnase.add(i)
        for i in range(int(crm.start), int(crm.stop)+1):
            len_regions += not exons[('dmel_'+crm.chrom, i)]
            num_snps += ('dmel_'+crm.chrom, i) in snps.index


    return num_snps, len_regions

def get_background_snprate(tsss, dnase, snps, exons):
    len_regions = 0
    num_snps = 0
    snprate_dict = {}
    for i in pb()(tsss.index):
        tss = tsss.ix[i]
        if isinstance(tss, pd.Series):
            n, l = get_snps_and_len(tss, dnase, snps, set(), exons)
            len_regions += l
            num_snps += n
            if l:
                snprate_dict[i] = n/l
        elif isinstance(tss, pd.DataFrame):
            seen_dnase2 = set()
            ns = 0
            ls = 0
            for i, tss in tss.iterrows():
                n, l = get_snps_and_len(tss, dnase, snps, seen_dnase2, exons)
                len_regions += l
                ls += l
                num_snps += n
                ns += n
            if ls:
                snprate_dict[i] = ns/ls
    return num_snps / len_regions, num_snps, len_regions, snprate_dict




if __name__ == "__main__":
    dnase = pd.read_table('Reference/BDTNP_R4_R6.tsv', na_values=['?'])
    tsss = pd.read_table('Reference/tss', index_col='fbgn')
    snps = pd.read_table('analysis/on_mel/melsim_variant.bed', names=['chrom', 'pos', '_', 'melsim'], header=None, index_col=[0,1])

    peak = pd.read_table('analysis/results/peak.tsv', index_col=0)
    logist = pd.read_table('analysis/results/logist.tsv', index_col=0)

    exons = defaultdict(bool)

    for i, entry in enumerate(open('Reference/mel_good.gtf')):
        entry = entry.split()
        chrom = entry[0]
        if entry[2] != 'exon': continue
        start = int(entry[3])
        stop = int(entry[4])
        for i in range(start, stop):
            exons[(chrom, i)] = True
        if i%1e4 == 0:
            print('.', end='')
            stdout.flush()
    print("Done loading exons")

    snp_rates = {}
    len_regions_dict = {}
    num_snps_dict = {}

    for gene in pb()(peak.index.union(logist.index)):
        tss = tsss.ix[gene]
        len_regions = 0
        num_snps = 0
        seen_dnase = set()
        if isinstance(tss, pd.Series):
            #print("Single CRM for gene", gene)
            num_snps, len_regions = get_snps_and_len(tss, dnase, snps, set(), exons)
        elif isinstance(tss, pd.DataFrame):
            #print("Multiple CRMs for gene", gene)
            for i, tss in tss.iterrows():
                n, l = get_snps_and_len(tss, dnase, snps, seen_dnase, exons)
                len_regions += l
                num_snps += n
        else:
            assert False
        if len_regions:
            len_regions_dict[gene] = len_regions
            num_snps_dict[gene] = num_snps
            snp_rates[gene] = num_snps/len_regions
        else:
            print('No DNase accessible regions near', gene)



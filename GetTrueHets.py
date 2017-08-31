from __future__ import print_function
import pandas as pd
from sys import argv
from collections import defaultdict, Counter
from multiprocessing import Pool
from progressbar import ProgressBar

def get_snp_counts(fname):
    snps_all = defaultdict(Counter)
    snps_i = pd.read_table(fname)
    for i, row in snps_i.iterrows():
        snps_current = snps_all[(row.CHR, row.POSITION)]
        for base, count in zip('ACGT', row['POS_A|C|G|T'].split('|')):
            snps_current[base] += int(count)
        # Note that the bases have been complemented here
        for base, count in zip('ACGT', row['NEG_A|C|G|T'].split('|')):
            snps_current[base] += int(count)
    return snps_all


if __name__ == "__main__":
    with Pool() as p:
        all_snp_files = [p.apply_async(get_snp_counts, (fname, ))
                         for fname in argv[1:]]


        snps_all = defaultdict(Counter)
        for snps_in_file in ProgressBar()(all_snp_files):
            snps_in_file = snps_in_file.get()
            for snp in snps_in_file:
                snps_all[snp] += snps_in_file[snp]

    print("Writing")
    with open('true_hets.tsv', 'w') as out:
        for key, counts in snps_all.items():
            most_common = counts.most_common(2)
            if len(most_common) < 2:
                continue
            _, c2 = most_common[1]
            if c2 > 10:
                print(*key, sep='\t', file=out)






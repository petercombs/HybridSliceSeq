from __future__ import print_function
import pandas as pd
from sys import argv
from collections import defaultdict, Counter

if __name__ == "__main__":
    snps_all = defaultdict(Counter)
    for fname in argv[1:]:
        snps_i = pd.read_table(fname)
        for i, row in snps_i.iterrows():
            snps_current = snps_all[(row.CHR, row.POSITION)]
            for base, count in zip('ACGT', row['POS_A|C|G|T'].split('|')):
                snps_current[base] += int(count)
            # Note that the bases have been complemented here
            for base, count in zip('TGCA', row['NEG_A|C|G|T'].split('|')):
                snps_current[base] += int(count)

    with open('true_hets.tsv', 'w') as out:
        for key, counts in snps_all.items():
            _, c2 = counts.most_common(2)[1]
            if c2 > 10:
                print(*key, sep='\t', file=out)






from __future__ import print_function
import pandas as pd
from argparse import ArgumentParser
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


def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--outfile', '-o', default='true_hets.tsv')
    parser.add_argument('--min-counts', '-m', type=int, default=20)
    parser.add_argument('snp_counts', nargs='*')
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    with Pool() as p:
        all_snp_files = [p.apply_async(get_snp_counts, (fname, ))
                         for fname in args.snp_counts]


        snps_all = defaultdict(Counter)
        for snps_in_file in ProgressBar()(all_snp_files):
            snps_in_file = snps_in_file.get()
            for snp in snps_in_file:
                snps_all[snp] += snps_in_file[snp]

    print("Writing")
    with open(args.outfile, 'w') as out:
        for key, counts in sorted(snps_all.items()):
            most_common = counts.most_common(2)
            if len(most_common) < 2:
                continue
            _, c2 = most_common[1]
            if c2 >= args.min_counts:
                print(*key, sep='\t', file=out)






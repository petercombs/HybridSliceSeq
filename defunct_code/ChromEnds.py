from __future__ import print_function
import pandas as pd
import Utils as ut
from sys import stderr

n = 50
sn = str(n)

if __name__ == "__main__":
    if 'chrom_of' not in locals():
        chrom_of = ut.get_chroms()
        coord_of = ut.get_coords()

        chrom_of_mel = (
            chrom_of
            .select(ut.not_contains(('\\', 'FBgn', ':')))
            .select(ut.not_startswith(('CR', 'mir-', 'mRpS', 'mRpL', 'RpS',
                                       'RpL')))
        )
        coord_of_mel = coord_of.ix[chrom_of_mel.index]
        ase = pd.read_table('analysis_godot/ase_summary_by_read.tsv',
                            **ut.pd_kwargs)
        maternals = {line.strip() for line in
                     open('analysis/results/strict_maternal_gene_names.txt')}
        ase_nomat = ase.drop(maternals)
        has_ase = ase_nomat.index[ase_nomat.count(axis=1) > 10]

    for chrom in ['2L', '2R', '3L', '3R', 'X']:
        genes_on_chrom = chrom_of_mel.index[chrom_of_mel ==
                                            chrom].intersection(has_ase)
        coords_on_chrom = coord_of_mel[genes_on_chrom].sort_values()
        print(chrom+'_left'+sn, 'POS',
              *coords_on_chrom.index[:n], sep='\t')
        print(chrom+'_right'+sn, 'POS',
              *coords_on_chrom.index[-n:], sep='\t')
        print(chrom+'_left_5mb', 'POS',
              *coords_on_chrom.index[coords_on_chrom < 5e6], sep='\t')
        print(chrom+'_right_5mb', 'POS',
              *coords_on_chrom.index[coords_on_chrom > coords_on_chrom.max() - 5e6], sep='\t')
        print(chrom+'_all', 'POS',
              *coords_on_chrom.index, sep='\t')
        print(chrom, len(coords_on_chrom), file=stderr)


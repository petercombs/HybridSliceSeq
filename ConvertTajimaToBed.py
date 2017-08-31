from __future__ import print_function
import pandas as pd
from sys import argv

tajima = pd.read_table(argv[1])

out = open(argv[2], 'w')

#out.write('#HEADER LINE\n')
for i, row in tajima.iterrows():
    print('chr{}\t{}\t{}\t{}'.format(
        row.CHROM,
        row.BIN_START, row.BIN_START+500,
        row.TajimaD),
        file=out,
        sep='\t',
    )


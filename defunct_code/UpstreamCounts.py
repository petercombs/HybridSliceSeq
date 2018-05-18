""" UpstreamCounts.py

Script to get the actual number of A's, C's, T's, and G's upstream of genes
"""
import pandas as pd
from Bio import SeqIO
from OrderedSeqRec import OrderedSeqRecord
from sys import argv
from collections import Counter
from progressbar import ProgressBar as pb


if __name__ == "__main__":
    coords = pd.read_table(argv[1])
    seqs = {rec.id: OrderedSeqRecord(rec) for rec in  SeqIO.parse(argv[2], 'fasta')}
    counts = Counter()
    for ix in pb()(coords.index):
        row = coords.ix[ix]
        counts.update(seqs[row.chrom][row.max_upstream:row.tss])

    print(counts)



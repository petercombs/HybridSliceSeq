#!/usr/bin/python
from __future__ import print_function
import pandas as pd
from argparse import ArgumentParser, FileType
from Bio import SeqIO
from Bio.Seq import MutableSeq

parser = ArgumentParser()
parser.add_argument('--emit-bed', '-b', default=None, type=FileType('w'))
parser.add_argument('fasta_file', type=open)
parser.add_argument('table', type=open)
parser.add_argument('outfile', type=str)
args = parser.parse_args()

seq_recs = [rec for rec in SeqIO.parse(args.fasta_file, 'fasta')]
sequence = {rec.id: MutableSeq(str(rec.seq)) for rec in seq_recs}
joint_vars = pd.read_table(args.table)

masked_sites = 0
corrected_sites = 0
for i, var in joint_vars.iterrows():
    if var['HOM-VAR'] == var['NCALLED'] and len(var['ALT']) == 1:
        sequence[var.CHROM][var.POS-1] = var['ALT']
        corrected_sites += 1
    elif var['HOM-REF'] != var['NCALLED']:
        if args.emit_bed and var['HOM-REF'] == 1 and var['HOM-VAR'] == 1:
            args.emit_bed.write('{}\t{}\t{}\t{}\n'.format(
                var.CHROM, 
                var.POS - 1,
                var.POS,
                var.ALT +"|"+sequence[var.CHROM][var.POS-1]
                )
                )
        for p in range(var.POS-1, var.POS + len(var.REF) - 1):
            sequence[var.CHROM][p] = 'N'
            masked_sites += 1
        

for rec in seq_recs:
    rec.seq = sequence[rec.id]

SeqIO.write(seq_recs, args.outfile, 'fasta')
print("Corrected sites: {:,}".format(corrected_sites))
print("Masked sites: {:,}".format(masked_sites))

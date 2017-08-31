#!/usr/bin/env python
from __future__ import print_function

from Bio import SeqIO
from sys import argv, exit
from os import path, makedirs

if len(argv) != 3:
    print("""Usage: {} FILENAME OUTDIR

FILENAME:       A FASTA file containing the sequences to split
OUTDIR:         The directory name to output split sequences to (will create if necessary)
""".format(argv[0]))
    exit(1)

if not path.exists(argv[2]):
    makedirs(argv[2])

if not path.isdir(argv[2]):
    print("{} is not a directory".format(argv[2]))
    exit(2) 

for rec in SeqIO.parse(argv[1], 'fasta'):
    rec.seq = rec.seq.upper()
    SeqIO.write(rec, path.join(argv[2], rec.id + '.fa'), 'fasta')


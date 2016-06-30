from Bio import SeqIO
from sys import argv

if __name__ == "__main__":

    fasta = {rec.id: rec for rec in SeqIO.parse(argv[1], 'fasta')}
    with open(argv[1].replace('fasta', 'pseq'), 'w') as outf:
        for recid in fasta:
            outf.write('{} \\{}\\\n'.format(recid, fasta[recid].seq))



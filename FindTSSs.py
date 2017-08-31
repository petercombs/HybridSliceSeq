#/usr/bin/env python
'''Quick script to find all TSSs in a GTF file


'''

from __future__ import print_function
from sys import stdin


def main():
    curr_transcript = None
    curr_chr = None
    curr_start = None
    curr_strand = None
    curr_gene = None
    curr_fbgn = None
    print('chr\tTSS_start\tgene_name\tfbgn\tfbtr')
    for line in stdin:
        data = line.strip().split('\t')
        chr = data[0]
        type = data[2]
        strand = data[6]
        min_coord = int(data[3])
        max_coord = int(data[4])
        annot = data[-1]
        annot = dict(item.replace('"', '').strip().split()
                     for item in annot.split(';')
                     if item.strip()
                    )

        if 'gene_name' not in annot: continue
        if type != 'exon': continue
        if annot['transcript_id'] != curr_transcript:
            if curr_strand == '-':
                print("{}\t{}\t{}\t{}\t{}"
                      .format(curr_chr, curr_start, curr_gene,
                              curr_fbgn, curr_transcript, ))
            curr_transcript = annot['transcript_id']
            curr_chr = chr
            curr_start = min_coord if strand == '+' else max_coord
            curr_strand = strand
            curr_gene = annot['gene_name']
            curr_fbgn = annot['gene_id']
            if strand == '+':
                print("{}\t{}\t{}\t{}\t{}"
                      .format(curr_chr, curr_start, curr_gene,
                             curr_fbgn, curr_transcript, ))






if __name__ == "__main__":
    main()


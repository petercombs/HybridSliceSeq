from __future__ import print_function
from sys import argv, stderr
from collections import defaultdict
from numpy import zeros
from progressbar import ProgressBar as pb
from argparse import ArgumentParser

# GTF Line:
#2L      FlyBase exon    7529    8116    .       +       .       transcript_id "FBtr0300689"; gene_id "FBgn0031208"; gene_name "CG11023";

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--max-distance', default=-1, type=int)
    parser.add_argument('--ordered', default=False, action='store_true')
    parser.add_argument('gtf_file')

    return parser.parse_args()



if __name__ == "__main__":
    args = parse_args()
    print(args, file=stderr)
    tss = defaultdict(dict)
    chrom_by_gene = {}
    gene_name_by_id = defaultdict(dict)
    strand_by_gene = {}
    in_exon = defaultdict(lambda : zeros(40000000, dtype=bool))
    last_on_chrom = {}
    print("Loading GFF", file=stderr)
    for line in open(args.gtf_file):
        chrom, _, ftype, low, hi, _, strand, _, desc = line.strip().split('\t')
        if ftype != 'exon': continue
        desc = {entr.strip().split(' ', 1)[0].strip() : 
                entr.strip().split(' ', 1)[1].strip('"')
                for entr in desc.split(';')
                if entr.strip()}
        gene_id = desc.get('gene_id', '???')
        gene_name = desc.get('gene_name', '???')
        transcript_id = desc.get('transcript_id', '???')
        if gene_id == '???':
            print( desc)
            raise ValueError("Misformed GFF Description: '{}'".format( line.strip().split('\t')))

        strand_by_gene[gene_id] = strand
        gene_name_by_id[chrom][gene_id] = gene_name
        chrom_by_gene[gene_id] = chrom
        low, hi = int(low), int(hi)
        if chrom not in last_on_chrom or hi > last_on_chrom[chrom][0]:
            last_on_chrom[chrom] = (hi, gene_id, gene_name)
        in_exon[chrom][low: hi+1] = True


        if strand == "+":
            if ((transcript_id not in tss[gene_id]) or (low < tss[gene_id][transcript_id])):
                tss[gene_id][transcript_id] = low
        elif strand == '-':
            if ((transcript_id not in tss[gene_id]) or (hi > tss[gene_id][transcript_id])):
                tss[gene_id][transcript_id] = hi
        else:
            assert False
    print("Finished Loading GFF", file=stderr)
    if args.ordered:
        print("chrom\tcoord_low\tcoord_hi\tFBgn\tgene_name\tstrand")
    else:
        print("chrom\ttss\tmax_upstream\tFBgn\tgene_name\tstrand")
    for gene in (gene
            for chrom in sorted(gene_name_by_id)
            for gene in sorted(gene_name_by_id[chrom], key=lambda x: min(tss[x].values()))
            ):
        strand = strand_by_gene[gene]
        chrom = chrom_by_gene[gene]
        name = gene_name_by_id[chrom][gene]
        if strand == '-' and last_on_chrom[chrom][1] == gene:
            gene_tss = sorted(tss[gene].values())
            prev_tss = gene_tss[1:] + ['end']
        else:
            if strand == '+':
                pos = min(tss[gene].values()) - 1
                while not in_exon[chrom][pos] and pos >= 0:
                    pos -= 1
                gene_tss = sorted(set(tss[gene].values()))
                prev_tss = [pos] + gene_tss[:-1]
            elif strand == '-':
                pos = max(tss[gene].values()) + 1
                while not in_exon[chrom][pos] and pos < len(in_exon[chrom]) - 1:
                    pos += 1
                gene_tss = sorted(set(tss[gene].values()))
                prev_tss = gene_tss[1:] + [pos]
        for curr_tss, prev_coord in zip(gene_tss, prev_tss):
            if args.max_distance > 0 and (prev_coord == 'end' or abs(curr_tss - prev_coord) > args.max_distance):
                if strand == '-':
                    prev_coord = curr_tss + args.max_distance
                elif strand == '+':
                    prev_coord = curr_tss - args.max_distance
                else:
                    raise ValueError("Bad strand value: " + strand)

            if args.ordered:
                curr_tss, prev_coord = sorted([curr_tss, prev_coord])
            print(
                    chrom,
                    curr_tss, prev_coord,
                    gene,
                    name,
                    strand,
                    sep='\t')




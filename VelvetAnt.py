"""Look for allele-specific splicing variants



Note that the name is a reference to the Mutlidae, which are wasps[1] that look
like ants (though not necessarily leaf cutter ants[2]).
[1] http://www.nature.com/nmeth/journal/v12/n11/full/nmeth.3582.html
[2] http://biorxiv.org/content/early/2016/03/16/044107
"""

from __future__ import print_function
from pysam import Samfile
from argparse import ArgumentParser, FileType
from collections import defaultdict
from os import path
from sys import stdout
from progressbar import ProgressBar
import pickle as pkl
import pandas as pd

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--snps-bed', '-S')
    parser.add_argument('--splicing-clusters', '-j',
                        help='''The output file from leafcutter_cluster''')
    parser.add_argument('--use-index', '-x', default=False, action='store_true')
    parser.add_argument('--outfile', '-o', default=stdout, type=FileType('w'))
    parser.add_argument('--min-reads', '-m', default=10,
                        help='Minimum number of reads to calculate a summary'
                        ' statistic')
    parser.add_argument('samfile', type=Samfile,
                       help='SAM/BAM alignment file to count')
    return parser.parse_args()

def get_phase(read, snps):
    phase = None
    for read_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
        if ref_pos + 1 in snps:
            if phase == None:
                try:
                    # 1 if alternate, -1 if reference
                    phase = -1 + 2*snps[ref_pos + 1].index(read.seq[read_pos])
                except ValueError:
                    return 0 # This SNP isn't in the dataset
            else:
                try:
                    new_phase = -1 + 2*snps[ref_pos + 1].index(read.seq[read_pos])
                except ValueError:
                    return 0
                if new_phase != phase:
                    return 0 # read seems misphased
    return phase



def get_snps(snpfile):
    snps = defaultdict(dict)
    if path.exists('true_hets.tsv'):
        true_hets = {tuple(line.strip().split()):True
                     for line in open('true_hets.tsv')
                    }
    else:
        true_hets = defaultdict(lambda x: True)
    if path.exists(snpfile+'.pkl'):
        return pkl.load(open(snpfile+'.pkl', 'rb'))

    if snpfile.endswith('.bed'):
        for line in open(snpfile):
            chrom, _, start, refalt = line.strip().split()
            if true_hets.get((chrom, start), True):
                snps[chrom][int(start)] = refalt.split('|')

    pkl.dump(snps, open(snpfile+'.pkl', 'wb'))
    return snps

def match_to_juncs(read, juncs, chrom=None, min_overlap=5):
    cigartuples = read.cigartuples
    # Note that in some cases, prefetching from the read is a lot faster than
    # always asking for read.cigartuples

    curpos = read.pos
    matched_juncs = []
    if chrom is None:
        chrom = read.reference_name
    for i, (cigar, length) in enumerate(cigartuples):
        if cigar == 3:
            # Skipped region from the reference
            if ((cigartuples[i-1][0] == 0 )
                and (cigartuples[i-1][1] >= min_overlap)
                and (cigartuples[i+1][0] == 0)
                and (cigartuples[i+1][1] >= min_overlap)):
                matched_juncs.append(juncs.get((chrom, curpos,
                                                 curpos+length+1)))
            curpos += length
        elif cigar == 0 or cigar == 1:
            # Match or insertion to reference
            curpos += length
    return matched_juncs

def load_juncs(juncfile):
    retval = defaultdict(dict)
    if juncfile.endswith('.gz'):
        import gzip
        opener = lambda x: gzip.open(x, 'rt')
    else:
        opener = open
    with opener(juncfile) as juncfile:
        next(juncfile)
        for line in juncfile:
            data = line.split()[0]
            chrom, left, right, clu = data.split(':')
            left = int(left)
            right = int(right)
            retval[clu][chrom, left, right] = data
    return retval


if __name__ == "__main__":
    args = parse_args()
    references = args.samfile.references
    juncs = load_juncs(args.splicing_clusters)
    snps = get_snps(args.snps_bed)
    all_juncs = {}
    for clu in juncs:
        all_juncs.update(juncs[clu])
    phases = {
        0: 'AMBIGUOUS',
        1: 'ALT',
        -1: 'REF',
        None: 'NO_SNPS'
    }
    out = pd.DataFrame(index=sorted(all_juncs.values())+[None], columns=phases.values(),
                       data=0)

    out.ix[None] = 0
    if args.use_index and args.samfile.has_index():
        for clu in ProgressBar()(juncs):
            left = 1e10
            right = 0
            for chrom, l, r in juncs[clu]:
                left = min(left, l)
                right = max(right, r)
            try:
                for read in args.samfile.fetch(chrom, left, right):
                    introns = match_to_juncs(read, juncs[clu])
                    if [i for i in introns if i]:
                        phase = phases[get_phase(read, snps[chrom])]
                        out.ix[introns, phase] += 1
            except ValueError:
                pass

    else:
        if args.samfile.has_index():
            num_reads = args.samfile.mapped
            pbar = ProgressBar(maxval=num_reads)
        else:
            pbar = lambda x: x

        for read in pbar(args.samfile):
            introns = match_to_juncs(read, all_juncs)
            #assert not (read.reference_name == 'dmel_2L' and 'N' in
                        #read.cigarstring)
            if [i for i in introns if i]:
                phase = phases[get_phase(read, snps[read.reference_name])]
                out.ix[introns, phase] += 1



    args.samfile.close()
    out['pref_index'] = (out.ALT - out.REF)/(out.ALT + out.REF)
    out.ix[out.ALT + out.REF < args.min_reads, 'pref_index'] = pd.np.nan
    out.index.name = 'clu'
    out.to_csv(args.outfile, sep='\t', na_rep='N/A')




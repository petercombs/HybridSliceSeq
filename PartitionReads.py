from os import path
from pysam import Samfile
from argparse import ArgumentParser
from collections import defaultdict
import pickle as pkl

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


def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--reference-species', '-r', default='ref')
    parser.add_argument('--alt-species', '-a', default='alt')
    parser.add_argument('--paired', '-p', default=False, action='store_true')
    parser.add_argument('snp_file')
    parser.add_argument('in_sam', type=Samfile)

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    if isinstance(args.in_sam.filename, str):
        fname = args.in_sam.filename
    elif isinstance(args.in_sam.filename, bytes):
        fname = args.in_sam.filename.decode()
    else:
        raise ValueError("Input filename is of type {}, which is not supported (yet)"
                         .format(type(fname)))
    base, ext = path.splitext(fname)
    print('{}_{}{}'.format(base, args.reference_species, ext),
                      'wb' if args.in_sam.is_bam else 'wr',
         )

    out_ref = Samfile('{}_{}{}'.format(base, args.reference_species, ext),
                      'wb' if args.in_sam.is_bam else 'wr',
                      template=args.in_sam)
    out_alt = Samfile('{}_{}{}'.format(base, args.alt_species, ext),
                      'wb' if args.in_sam.is_bam else 'wr',
                      template=args.in_sam)
    out_amb = Samfile('{}_{}{}'.format(base, 'ambig', ext),
                      'wb' if args.in_sam.is_bam else 'wr',
                      template=args.in_sam)

    snps = get_snps(args.snp_file)
    outs = [out_amb, out_alt, out_ref]

    unmatched_reads = [{}, {}]
    chroms = args.in_sam.references
    try:
        from progressbar import ProgressBar
        iterator = ProgressBar(maxval=args.in_sam.mapped)(args.in_sam,)
        print("Using progress bar")

    except Exception as err:
        print("We had an issue: {}".format(err))
        iterator = args.in_sam

    for read in iterator:
        if not args.paired:
            phase = get_phase(read, snps[chroms[read.reference_id]])
            if phase is None:
                phase = 0
            outs[phase].write(read)
        elif read.qname in unmatched_reads[read.is_read1]:
            other_read = unmatched_reads[read.is_read1].pop(read.qname)
            phase_self = get_phase(read, snps[chroms[read.reference_id]])
            phase_other = get_phase(other_read, snps[chroms[read.reference_id]])
            if read.is_read2:
                read, other_read = other_read, read
                phase_self, phase_other = phase_other, phase_self
            else:
                pass
            if ((phase_self == phase_other) and phase_self is not None) or (phase_self is not None and phase_other is None):
                outs[phase_self].write(read)
                outs[phase_self].write(other_read)
            elif phase_self is None and phase_other is not None:
                outs[phase_other].write(read)
                outs[phase_other].write(other_read)
            else:
                outs[0].write(read)
                outs[0].write(other_read)
        else:
            unmatched_reads[read.is_read2][read.qname] = read

    print("There are {} unmatched read 1s \n and {} unmatched read 2s"
          .format(len(unmatched_reads[0]), len(unmatched_reads[1])))

    for read in unmatched_reads[0].values():
        phase = get_phase(read, snps[chroms[read.reference_id]]) or 0
        outs[phase].write(read)

    for read in unmatched_reads[1].values():
        phase = get_phase(read, snps[chroms[read.reference_id]]) or 0
        outs[phase].write(read)


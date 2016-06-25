from __future__ import print_function
from os import path
from GetMapStats import count_reads

def parse_args():
    from argparse import ArgumentParser

    parser = ArgumentParser(description=("Take a collection of fpkm_tracking "
                                         "files and makes a summary table "
                                         "with all of the data as a TSV file"))
    parser.add_argument('--params', '-p', default=False,
                        dest='has_params',
                        help="Parameters file including renaming conventions: "
                        'Directory in column "Label", stage in column "Stage"')
    parser.add_argument('mapped_bamfile',
                        help='The bam file to look in for mapped reads')

    args = parser.parse_args()
    return args



if __name__ == "__main__":
    args = parse_args()
    args
    read_counts = count_reads(args.mapped_bamfile, has_carrier=False,
                              count_unique=True, count_all=True,
                              force_count=True)

    out = open(path.splitext(args.mapped_bamfile)[0]+'.mapstats', 'w')
    print("UniqueMapped\tAll_mapped\tTotal", file=out)
    print("{}\t{}\t{}".format(*read_counts), file=out)
    out.close()

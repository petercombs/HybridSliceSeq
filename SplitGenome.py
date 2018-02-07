from __future__ import print_function
from argparse import ArgumentParser
from os import path, makedirs

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--size', '-s', default=10e6, type=float )
    parser.add_argument('chrom_size_file')
    parser.add_argument('out_base')
    args = parser.parse_args()
    args.size = int(args.size)
    return args

if __name__ == "__main__":
    args = parse_args()
    i = args.size
    out_filenum = 0
    try:
        makedirs(path.dirname(args.out_base))
    except FileExistsError:
        pass
    out_file = open('{}_{}.bed'.format(args.out_base, out_filenum), 'w')
    for line in open(args.chrom_size_file):
        chr_start = 0
        chrom, chr_size = line.split()
        chr_size = int(chr_size)

        print(chrom, chr_size)
        while chr_size > chr_start + i:
            print(chrom, chr_start, chr_start + i, sep='\t', file=out_file)
            out_filenum += 1
            out_file.close()
            out_file = open('{}_{}.bed'.format(args.out_base, out_filenum), 'w')
            chr_start += i
            i = args.size
        print(chrom, chr_start, chr_size, sep='\t', file=out_file)
        i -= (chr_size - chr_start)




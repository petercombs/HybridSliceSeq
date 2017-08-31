from __future__ import print_function
import pysam
import gzip
from sys import argv
from os import path
import shutil
import tempfile
from argparse import ArgumentParser
import re

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--single-file-per-run', '-s', default=False,
                        action='store_true',
                        help='This flag indicates that the sequencing core '
                        'split files by run only, not within each run')
    parser.add_argument('--organism', default='D. melanogaster')
    parser.add_argument('--insert-size', default=-1)
    parser.add_argument('--insert-std', default=-1)
    parser.add_argument('bamfile')
    parser.add_argument('outbase')

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    sf = pysam.Samfile(args.bamfile)
    rfs = [entry for entry in
           sf.header['PG'][0]['CL'].split()
           if entry.endswith('.gz') or entry.endswith('.fastq')][0]
    rfs = sorted(rfs.split(','))

    cuff_data =  path.join(
        path.dirname(args.bamfile),
        'genes.e.log'
    )
    if (args.insert_size < 0 or args.insert_std < 0) and path.exists(cuff_data):
        cuff_data = re.findall('Estimated [^:]*:[^\n]*\n',
                               open(cuff_data).read())
        for line in cuff_data:
            if line.strip().startswith('Estimated Mean'):
                args.insert_size = line.split(':')[1].strip()
            if line.strip().startswith('Estimated Std Dev'):
                args.insert_std = line.split(':')[1].strip()
            if not (args.insert_size == -1 or args.insert_std == -1):
                break
        else:
            print("Can't find frag data for ", args.outbase)
    else:
        print("Can't find frag data for ", args.outbase)

    fnames = []
    if args.single_file_per_run:
        for i, rf in enumerate(rfs):
            extension = '.fastq' + rf.split('.fastq', 1)[1]
            outfname = '{base}_R1_run{i}{ext}'.format(
                base=args.outbase,
                i=i+1,
                ext=extension)

            shutil.copyfile(rf, outfname)
            fnames.append(outfname)
            rf_handle = gzip.open(rf, 'rt') if rf.endswith('.gz') else open(rf)
            next(rf_handle)
            readlen = len(next(rf_handle).strip())
            rf2 = rf.replace('R1', 'R2')
            outfname2 = outfname.replace('R1', 'R2')
            print('RF:',
                  outfname, 'fastq', '',
                  'Illumina NextSeq 500' if readlen< 80 else 'Illumina HiSeq 2500',
                  readlen,
                  'paired-end' if path.exists(rf2) else 'single',
                  sep='\t',
                 )

            if path.exists(rf2):
                print('RF:',
                      outfname2, 'fastq', '',
                      'Illumina NextSeq 500' if readlen< 80 else 'Illumina HiSeq 2500',
                      readlen,
                      'paired-end' if path.exists(rf2) else 'single',
                      sep='\t',
                     )
                print("PF:",
                      outfname, outfname2,args.insert_size, args.insert_std,
                     sep='\t')
                shutil.copyfile( rf2, outfname2)
                fnames.append(outfname2)
    else:
        raise NotImplementedError("I don't have a good way to do this yet")
        with open(argv[2]+'_R1.fastq.gz', 'wb') as outf:
            for i, rf in enumerate(rfs):
                if rf.endswith('.gz'):
                    with open(rf, 'rb') as readfile:
                        shutil.copyfileobj(readfile, outf)
                else:
                    with gzip.open(tempfile.TemporaryFile(), 'w+b') as reads:
                        for line in open(rf):
                            reads.write(line)
                        reads.seek(0)
                        shutil.copyfileobj(reads, outf)

        if any([path.exists(rf.replace('R1', 'R2')) for rf in rfs]):
            outf2 = gzip.open(argv[2]+'R2.fastq.gz', 'w')
            for rf in rfs:
                rf = rf.replace('R1', 'R2')
                if not path.exists(rf): continue

                if rf.endswith('.gz'):
                    infile = gzip.open(rf)
                else:
                    infile = open(rf)
                for line in infile:
                    outf2.write(line)
            outf2.close()


    print("SN:",
          "Sample N", args.outbase, args.outbase.split('_sl')[0],
          args.organism,
          'embryonic blastoderm'
          'cycle 14C',
          'mRNA',
          '',
          '',
          *fnames,
          sep='\t')



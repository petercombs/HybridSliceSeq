from __future__ import division, print_function
import pandas as pd
import gzip
import multiprocessing as mp
from itertools import repeat
from os import path
from glob import glob
from pysam import Samfile
from collections import defaultdict

def parse_args():
    from argparse import ArgumentParser

    parser = ArgumentParser(description=("Take a collection of fpkm_tracking "
                                         "files and makes a summary table "
                                         "with all of the data as a TSV file"))
    parser.add_argument('--confidence', '-c', default=False,
                        dest='conf', action='store_true',
                        help="Include confidence intervals")
    parser.add_argument('--params', '-p', default=False,
                        dest='has_params',
                        help="Parameters file including renaming conventions: "
                        'Directory in column "Label", stage in column "Stage"')
    parser.add_argument('--key', '-k', default='gene_short_name',
                        help='The column to combine on (FBgn in tracking_id)')
    parser.add_argument('--mapped-bamfile', '-b', default='assigned_dmelR.bam',
                        help='The bam file to look in for mapped reads')
    parser.add_argument('--in-subdirectory', default=None,
                        help='Subdirectory in '
                        'basedir/sample/subdirectory/genes.fpkm_tracking')
    parser.add_argument('--out-basename', '-o', dest='basefile',
                        default='summary',
                        help='The base of the output filename to which'
                        ' modifiers may be appended depending on options'
                        '(defaults to "summary")')
    parser.add_argument('--count-unique', '-u', default=False,
                        action='store_true',
                        help='When removing samples, use the number of unique '
                        'reads, not total number of mappings')
    parser.add_argument('--count-all', '-a', default=False,
                        action='store_true',
                        help='Count all reads in FASTQ files.  This will be '
                        'more accurate, but also potentially slower than using '
                        'the heuristic that all but the last files are capped at '
                        '4M reads')
    parser.add_argument('--translate-labels', '-T', default=False,
                        action='store_true')
    parser.add_argument('basedir',
                        help='The directory containing directories, which '
                        'contain genes.fpkm_tracking files')

    args = parser.parse_args()
    return args


def get_stagenum(name, series, dir):
    if not [c for c in name if c.isdigit()]:
        return 0
    # Slice until the first digit
    name_base = name[:[i for i, c in enumerate(name) if c.isdigit()][0]]
    dir = {'+': 1, '?': 1, '-': -1}[dir]
    return (sorted({ix for ix in series if ix.startswith(name_base)})[::dir]
            .index(name)) + 1

def mp_count_reads(args):
    return count_reads(*args)

def count_reads(bamname, has_carrier=False, count_unique=False, count_all=False,
               force_count=False):
    mapfile = path.splitext(bamname)[0]+'.mapstats'
    if path.exists(mapfile) and not force_count:
        return list(pd.read_table(mapfile).ix[0])
    sf = Samfile(bamname)
    u_mapped = pd.np.nan
    if count_unique:
        u_mapped = 0
        for read in sf:
            u_mapped += not read.is_secondary

    reads = 0
    rfs = [entry for entry in
           sf.header['PG'][0]['CL'].split()
           if entry.endswith('.gz') or entry.endswith('.fastq')][0]
    rfs = rfs.split(',')
    rfs_by_dir = defaultdict(list)
    for rf in rfs:
        rfs_by_dir[path.dirname(rf)].append(rf)

    reads = pd.np.nan
    if not has_carrier:
        reads = 0
        for rfs in rfs_by_dir.values():
            rfs = sorted(rfs)
            if count_all:
                for rf in rfs:
                    for i, line in enumerate(gzip.open(rfs[-1])):
                        pass
                    reads += (i+1)//4
            else:
                reads += 4e6 * (len(rfs) - 1)
                for i, line in enumerate(gzip.open(rfs[-1])):
                    pass
                reads += (i+1)//4
    return u_mapped, sf.mapped, reads

if __name__ == "__main__":
    args = parse_args()
    if args.in_subdirectory:
        fnames = glob(path.join(args.basedir, '*', args.in_subdirectory,
                                args.mapped_bamfile))
    else:
        fnames = glob(path.join(args.basedir, '*', args.mapped_bamfile))
    has_carrier = defaultdict(bool)
    if args.has_params:
        params = (pd.read_table(args.has_params,
                                    comment='#',
                                    converters={'Label': str},
                                    #keep_default_na=False, na_values='---',
                                   )
                  )
        params = params.dropna(how='any')
        for dirname in params.Label:
            has_carrier[dirname] = (set(params.ix[params.Label == dirname,
                                                  'CarrierSpecies'])
                                    != set(['---']))
            has_carrier[dirname] = False
    else:
        has_carrier = defaultdict(bool)

    pool = mp.Pool(30)
    fdirs = [fname.split('/')[1] for fname in fnames]
    res = pool.map(mp_count_reads,
                   zip(fnames,
                       [has_carrier[dirname] for dirname in fdirs],
                       repeat(args.count_unique),
                       repeat(args.count_all)))

    if args.has_params:
        new_fdirs = []
        for dirname in fdirs:
            indices = params.Label.index[params.Label == dirname]
            if len(indices):
                i = indices[0]
                new_dirname = "{genotype}_cyc{stage}_sl{num:02}".format(
                    genotype=params.ix[i, 'SampleGenotype'],
                    stage=params.ix[i, 'Stage'],
                    num=get_stagenum(dirname, params.Label,
                                     params.ix[i, 'Direction']))
            else:
                new_dirname = dirname
            new_fdirs.append(new_dirname)
    df = pd.DataFrame(index=fdirs,
                      data=pd.np.array(res),
                      columns=['UniqueMapped', 'AllMapped', 'RawReads'],
                     )
    if args.translate_labels:
        df['LongName'] = new_fdirs
    df.sort_index().to_csv(path.join(args.basedir,
                                     'map_stats.tsv'),
                           sep='\t',
                          )



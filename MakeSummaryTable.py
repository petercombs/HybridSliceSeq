"""MakeSummaryTable

Takes a collection of fpkm_tracking files and makes a summary table with all of
the data as a CSV file.  Arguments:
    1) superdirectory
    2) -c -- also include confidence intervals

"""
from __future__ import division, print_function
import pandas
from os import path
from glob import glob
from pysam import Samfile
import gzip
import sys


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
    parser.add_argument('--strip-low-reads', '-s', default=0, type=int,
                        help='Remove samples with fewer than N counts (off by'
                        'default)')
    parser.add_argument('--strip-on-unique', '-u', default=False,
                        action='store_true',
                        help='When removing samples, use the number of unique '
                        'reads, not total number of mappings')
    parser.add_argument('--strip-as-nan', '-n', default=False,
                        action='store_true',
                        help='When stripping a sample, replace all data with'
                        ' NaN')
    parser.add_argument('--strip-low-map-rate', '-m', default=0, type=float,
                        help='Remove samples with less than X%% of reads '
                        "mapping (off by default)")
    parser.add_argument('--mapped-bamfile', '-b', default='assigned_dmelR.bam',
                        help='The bam file to look in for mapped reads')
    parser.add_argument('--in-subdirectory', default=None,
                        help='Subdirectory in '
                        'basedir/sample/subdirectory/genes.fpkm_tracking')
    parser.add_argument('--filename', default='genes.fpkm_tracking',
                        help=('Filename of the per-sample gene '
                              'expression table'))
    parser.add_argument('--column', '-C', default='FPKM',
                        help='Column to read out (either name or number)')
    parser.add_argument('--no-header', dest='header', action='store_false',
                        default=True,
                        help='No header line in the file')
    parser.add_argument('--out-basename', '-o', dest='basefile',
                        default='summary',
                        help='The base of the output filename to which'
                        ' modifiers may be appended depending on options'
                        '(defaults to "summary")')
    parser.add_argument('--map-stats', default=None,
                        help="File containing mapping statistics.")
    parser.add_argument('--float-format', default='%8.2f',
            help="Floating point format string")
    parser.add_argument('--make-geo', default=False, action='store_true',
            help="Create a makefile to easily package files for GEO submission")
    parser.add_argument('--skip-bad', default=False, action='store_true',
                        help="Try to skip files with no readable data")
    parser.add_argument('basedir',
                        help='The directory containing directories, which '
                        'contain genes.fpkm_tracking files')

    args = parser.parse_args()
    try:
        args.column = int(args.column)
    except ValueError:
        pass
    try:
        args.key = int(args.key)
    except ValueError:
        pass
    if args.map_stats:
        try:
            args.map_stats = pandas.read_table(args.map_stats, index_col=0)
        except:
            args.map_stats = None

    if args.strip_low_map_rate:
        args.strip_on_unique = True
        args.strip_low_reads = max(args.strip_low_reads, 1)
    return args


from GetMapStats import get_stagenum

def get_expr_values(fname):
    try:
        table = pandas.read_table(fname, na_values=['-', 'NA', 'N/A'], converters={args.key: str},
                                  comment='#', skip_blank_lines=True,
                              keep_default_na=False,
                              header=None if not args.header else 0)
    except Exception as err:
        print("Error on", fname, file=sys.stderr)
        print("Error on", fname, file=sys.stdout)
        sys.stderr.flush()
        sys.stdout.flush()
        if args.skip_bad:
            return (None, None)
        else:
            raise err
    alldir, fname = path.split(fname)
    if args.in_subdirectory:
        alldir = (alldir.replace(args.in_subdirectory, '')
                  .replace('//', '/')
                  .strip('/'))
    basedir, dirname = path.split(alldir)
    old_dirname = dirname
    if args.column not in table.columns:
        print("Couldn't find {} in table {}".format(
            args.column, path.join(alldir, fname)
        ))
        return (None, None)
        raise ValueError("argument {} not in table {} columns: {}".format(
            args.column, path.join(alldir, fname), table.columns))
    if args.key in table.columns:
        table = (table.drop_duplicates(args.key)
                 .drop(labels=['ORIENTATION'], axis=1, errors='ignore')
                 .dropna(axis=1, how='all')
                 .dropna(axis=0, how='any'))
    if args.key not in table.index and args.key not in table.columns:
        print("Error on", old_dirname, fname, file=sys.stderr)
        print("Error on", old_dirname, fname, file=sys.stdout)
        print(table, file=sys.stderr)
        sys.stderr.flush()
        sys.stdout.flush()
        return (None, None)
    if isinstance(args.key, int) and args.key not in table.columns:
        table.set_index(table.columns[args.key], inplace=True,
                        verify_integrity=True)
    else:
        table.set_index(args.key, inplace=True, verify_integrity=True)
    if args.has_params and dirname not in params.index:
        return (None, None)

    if args.has_params:
        new_dirname = "{genotype}_cyc{stage}_sl{num:02}".format(
            genotype=params.ix[dirname, 'SampleGenotype'],
            stage=params.ix[dirname]['Stage'],
            num=get_stagenum(dirname, params.index,
                             params.ix[dirname, 'Direction']))
        print(dirname, '=', new_dirname)
        dirname = new_dirname
    else:
        print(dirname)

    if args.column not in table.columns:
        return (dirname, pandas.np.nan)
    skip = False
    if args.strip_low_reads:
        if (args.map_stats is not None) and dirname in args.map_stats.index:
            col = 'UniqueMapped' if args.strip_on_unique else 'AllMapped'
            reads = args.map_stats.ix[dirname, col]
        elif (args.map_stats is not None) and old_dirname in args.map_stats.index:
            col = 'UniqueMapped' if args.strip_on_unique else 'AllMapped'
            reads = args.map_stats.ix[old_dirname, col]
        else:
            if (args.map_stats is not None):
                print("Missing {} in mapping stats".format(dirname))
            try:
                sf = Samfile(path.join(alldir, args.mapped_bamfile))
                if args.strip_on_unique:
                    reads = 0
                    for read in sf:
                        reads += not read.is_secondary
                        if ((reads > args.strip_low_reads)
                            and not args.strip_low_map_rate):
                            break
                else:
                    reads = sf.mapped
            except IOError:
                print("Error reading", path.join(alldir, args.mapped_bamfile))
                reads = 0
        skip += reads < args.strip_low_reads
    if (args.strip_low_map_rate
        and not skip
        and not (args.has_params
                 and 'CarrierSpecies' in params.columns
                 and params.ix[old_dirname, 'CarrierSpecies'] != '---')
       ):

        if (not (args.map_stats is None)
            and (dirname in args.map_stats.index)
            and (args.map_stats.ix[dirname, 'RawReads'] > 0)):
            total_reads = args.map_stats.ix[dirname, 'RawReads']
        elif (not (args.map_stats is None)
              and (old_dirname in args.map_stats.index)
              and (args.map_stats.ix[old_dirname, 'RawReads'] > 0)):
            total_reads = args.map_stats.ix[old_dirname, 'RawReads']
        else:
            rfs = [entry for entry in
                   sf.header['PG'][0]['CL'].split()
                   if entry.endswith('.gz') or entry.endswith('.fastq')][0]
            rfs = sorted(rfs.split(','))
            total_reads = 4e6 * (len(rfs) - 1)
            for i, line in enumerate(gzip.open(rfs[-1])):
                pass
            total_reads += i//4
        map_rate = reads / total_reads
        if map_rate > 1: # Quick and dirty correction for paired end reads
            map_rate /= 2
        skip += map_rate < (args.strip_low_map_rate / 100)
        if skip:
            print(reads, total_reads, map_rate,
                  args.strip_low_map_rate / 100)

    if skip:
        if args.strip_as_nan:
            from numpy import nan
            print("NaNing", dirname, "\t{:,} reads".format(reads))
            table.ix[:] = nan
        else:
            print("Skipping", dirname)
            return None, None
    if args.conf:
        return (
            dirname+"_"+args.column,
            table.ix[:, args.column],
            dirname+"_"+args.column+"_conf_lo",
            table.ix[:, args.column+"_conf_lo"],
            dirname+"_"+args.column+"_conf_hi",
            table.ix[:, args.column+"_conf_hi"],
               )
    else:
        try:
            ret = table.ix[:, args.column]
        except Exception as e:
            print("Problem with ", old_dirname)
            raise e
        return (dirname+"_"+ret.name,ret )

if __name__ == "__main__":
    args = parse_args()
    if args.in_subdirectory:
        fnames = glob(path.join(args.basedir, '*', args.in_subdirectory,
                                args.filename))
    else:
        fnames = glob(path.join(args.basedir, '*', args.filename))
    if args.has_params:
        params = (pandas.read_table(args.has_params,
                                    comment='#',
                                    converters={'Label': str},
                                    #keep_default_na=False, na_values='---',
                                   )
                  .drop_duplicates(subset=['Label']))
        params.set_index('Label', inplace=True)
        params = params.dropna(how='any')


    import multiprocessing as mp
    pool = mp.Pool(30)
    res = pool.map(get_expr_values, fnames)
    if args.conf:
        names, cols, names_lo, cols_lo, names_hi, cols_hi = zip(*res)
        df = pandas.DataFrame(dict(zip(names, cols))
                              +dict(zip(names_lo, cols_lo))
                              +dict(zip(names_hi, cols_hi)))
    else:
        names, cols = zip(*res)
        df = pandas.DataFrame(dict(zip(names, cols)))


    df.sort_index(axis=1).to_csv(path.join(args.basedir,
                                           args.basefile
                                           + ('_in_{}'.format(args.in_subdirectory)
                                              * bool(args.in_subdirectory))
                                           + ('_with_conf' * args.conf)
                                           + '.tsv'),
                                 float_format=args.float_format,
                                 sep='\t', na_rep='---')

    if args.make_geo:
        if not path.exists(path.join(args.basedir, 'geo')):
            import os
            os.makedirs(path.join(args.basedir, 'geo'))

        commands = open(path.join(args.basedir, 'geo.make'), 'w')
        print("Names: '{}'".format(names))
        commands.write('all: {} \n'.format(' '.join('geo/'+
                                                    f.replace('_FPKM',
                                                               '_R1.fastq.gz.md5')
                                                     for f in names if f)))
        commands.write('%.md5 : % \n\tmd5sum $< > $@\n\n\n')
        for bamname, fname  in zip(fnames, names):
            try:
                bamname = bamname.replace(args.filename, args.mapped_bamfile)
                commands.write(('geo/{fname}_R1.fastq.gz: {bamname}\n'
                    '\tpython CompileForGEO.py {bamname} geo/{fname}\n\n')
                    .format(fname=fname.replace('_FPKM', ''), bamname=bamname))
            except Exception as err:
                print("Skipping exception: \n", bamname, fname, err)

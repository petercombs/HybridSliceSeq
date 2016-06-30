#!/usr/bin/env python
from __future__ import division, print_function
from collections import defaultdict
import pysam
import re
import sys
from glob import glob
from os import path
from scipy import stats
from numpy import array, log, exp
import progressbar as pbar

if len(sys.argv) < 3:
    print("""Usage: python CheckCoverage.py <GTF-File> BAMfile [BAMfile ...]""")
    sys.exit(1)

gtf_fname = sys.argv[1]  # 'Reference/dmel-all-r5.42.gtf'
analysis_dir = 'analysis'

starts = set()
curr_len = -1
coverage = 0

cutoff = 0


def analyze_bamfile(bam_fname):
    print(bam_fname,)
    bam_file = pysam.Samfile(bam_fname, 'rb')

    coverages = defaultdict(lambda: [0, 0, set()])
    parent = ''
    curr_len = 0
    coverage = 0
    starts = set()

    f = open(gtf_fname)
    for i, line in enumerate(f):
        pass
    pb = pbar.ProgressBar(widgets=[bam_fname, pbar.Bar(), pbar.ETA()],
                          maxval=i+1)
    f = open(gtf_fname)
    for line in pb(f):
        if line.startswith('#'):
            continue
        if line.startswith('>'):
            break
        data = line.split()
        chrom = data[0]
        kind = data[2]
        start = int(data[3]) - 1
        stop = int(data[4])
        fbtr_finder = re.compile('FBtr[0-9]*')
        # parent = fbtr_finder.findall(line)[0]

        if kind == 'exon':
            fbtrs = fbtr_finder.findall(line)
            if not fbtrs:
                continue
            fbtr = fbtrs[0]
            coverages[fbtr][1] += (stop - start)

            starts = set()
            coverage = 0
            for read in bam_file.fetch(chrom, start, stop):
                coverage += 1
                starts.add(read.pos)
            coverages[fbtr][0] += coverage
            coverages[fbtr][2].update(starts)

    pb.finish()
    curr_lens, rpks, uniques = zip(*coverages.values())
    dir, fname = path.split(bam_fname)

    return dir, rpks, uniques, curr_lens

if __name__ == "__main__":
    import multiprocessing as mp

    POOL = mp.Pool(20)
    res = POOL.map(analyze_bamfile,
                   [f for f in sys.argv[2:] if f.endswith('.bam')])
    all_dirs, all_rpks, all_pct_uniques, all_lens = zip(*res)

    import pickle
    out_fh = open('checkcoverage.pkl', 'wb')
    pickle.dump({'dirs': all_dirs, 'rpks': all_rpks,
                 'pct_uniques': all_pct_uniques, 'lens': all_lens},
                out_fh)
    for fname, rpks, uniques, curr_lens in res:
        print(fname)
        try:
            xs = array(rpks)
            ys = array([len(u)/(curr_len + 1)
                       for u, curr_len in zip(uniques, curr_lens)])
            cutoff = max(xs[ys < .1])
            reg = stats.linregress(log(xs[(xs < cutoff)*(xs > 0)*(ys > 0)]),
                                   log(ys[(xs < cutoff)*(xs > 0)*(ys > 0)]))
            print("exp(%f) * x ** %f" % (reg[1], reg[0]))
            print("Duplicate badness score: ", exp(-reg[1]-.38))
        except Exception as exc:
            print(exc)

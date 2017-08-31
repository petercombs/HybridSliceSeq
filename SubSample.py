from __future__ import print_function
from pysam import Samfile
from numpy.random import random
from heapq import heappush, heappushpop
from os import path, makedirs
from sys import argv, stdout


def main():

    fname = argv[1]
    millions = list(map(lambda x: 1e6 * float(x), argv[2:]))
    print(millions)
    stdout.flush()
    return subsample(fname, millions, paired=True)


def subsample(fn, ns=None, paired=False):
    if ns is None:
        fn, ns = fn
    sample = []
    count = 0
    outdir_base = path.join(path.dirname(fn), 'subset')
    sf = Samfile(fn)
    try:
        i_weight = float(sf.mapped)/max(ns)
        print("Read out ", i_weight)
    except ValueError:
        i_weight = 0.0
        for read in sf:
            i_weight += 1
        print("Counted ", i_weight)
        i_weight /= float(max(ns))
        sf = Samfile(fn)

    if paired:
        read_2s = {}
    print(fn, count, i_weight)
    for i, read in enumerate(sf):
        key = random()**i_weight
        if not paired or read.is_read1:
            if len(sample) < max(ns):
                heappush(sample, (key, i+count, read))
            else:
                dropped = heappushpop(sample, (key, i+count, read))
                if paired:
                    read_2s.pop(dropped[-1].qname, None)
        elif paired:
            read_2s[read.qname] = read
        else:
            assert ValueError("I don't know how we got here")


    count += i

    for n in ns:
        outdir = outdir_base + '{:04.1f}M'.format(n/1e6)
        try:
            makedirs(outdir)
        except OSError:
            pass
        sampN = sorted(sample, reverse=True)[:int(n)]
        print("Kept {: >12,} of {: >12,} reads".format(len(sampN), count))
        print(fn, '->', outdir)
        stdout.flush()
        of = Samfile(path.join(outdir, 'accepted_hits.bam'),
                     mode='wb', template=sf)
        sample.sort(key=lambda heap_item: (heap_item[-1].tid, heap_item[-1].pos))
        missing_mates = 0
        for key, pos, read in sampN:
            of.write(read)
            if paired and read.is_proper_pair:
                if read.qname not in read_2s:
                    missing_mates += 1
                    continue
                of.write(read_2s[read.qname])
        of.close()
    sf.close()
    print(missing_mates)
    return [count for key, read, count in sample]

if __name__ == "__main__":
    subset_results = main()

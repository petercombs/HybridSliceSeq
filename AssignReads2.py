import pysam
import sys
from os import path
from collections import defaultdict, Counter
from progressbar import ProgressBar, ETA, Bar, Percentage
from argparse import Namespace


class my_defaultdict(dict):
    def __init__(self, default_factory, basename, other_args):
        self.default_factory = default_factory
        self.basename = basename
        self.other_args = other_args

    def __missing__(self, key):
        self[key] = value = self.default_factory(self.basename % key,
                                                 **self.other_args)
        return value


def get(read, tag):
    return {tname.upper(): val for tname, val in read.tags}[tag.upper()]


def get_nh(read):
    return {tag.upper(): val for tag, val in read.tags}['NH']


def get_species(read):
    return references[read.rname].split('_')[0]


def process_read(read):
    if not read.tags:
        print read
        print read.tags
        assert False
    nh = get_nh(read)
    species = get_species(read)
    if nh == 1:
        # assigned.write(read)
        specific_files[species].write(read)
        species_counts[species] += 1
        return
    else:
        resolve_multiread(read, nh, species)


def resolve_multiread(read, nh, species):
    nm = get(read, 'NM')
    has_multi_frags = bool(0x1 & read.flag)
    has_multi_frags = True
    if has_multi_frags and read.is_read2:
        to_be_resolved_counts2[read.qname] += 1
        if nm < to_be_resolved_vals2[read.qname][species]:
            to_be_resolved_vals2[read.qname][species] = nm
            to_be_resolved_reads2[read.qname][species] = read

        if nh == to_be_resolved_counts2[read.qname]:
            dbs = Namespace()
            dbs.to_be_resolved_vals = to_be_resolved_vals2
            dbs.to_be_resolved_counts = to_be_resolved_counts2
            dbs.to_be_resolved_reads = to_be_resolved_reads2
            on_last_multiread(dbs, read)

    elif has_multi_frags and read.is_read1:
        to_be_resolved_counts[read.qname] += 1
        if nm < to_be_resolved_vals[read.qname][species]:
            to_be_resolved_vals[read.qname][species] = nm
            to_be_resolved_reads[read.qname][species] = read

        if nh == to_be_resolved_counts[read.qname]:
            dbs = Namespace()
            dbs.to_be_resolved_vals = to_be_resolved_vals
            dbs.to_be_resolved_counts = to_be_resolved_counts
            dbs.to_be_resolved_reads = to_be_resolved_reads
            on_last_multiread(dbs, read)
    else:
        pass
        # print "didi we not have multiple frags?"
        # print has_multi_frags
        # print read.is_read1, read.is_read2
        # print read.qname
        # assert False
        #  WTF are we doign here?


def on_last_multiread(dbs, read):
    # Sort out the reads
    if len(dbs.to_be_resolved_reads[read.qname]) == 1:
        # Looks like there were multiple hits from the same species.
        # Report the best, or if equal quality, the first (which
        # tophat would've given anyways)
        species = get_species(read)
        # assigned.write(read)
        specific_files[species].write(read)
        species_counts[species] += 1
    else:
        # Hits from multiple species
        vals = sorted([(val, spec) for spec, val in
                       dbs.to_be_resolved_vals[read.qname].iteritems()])
        diff_val = vals[1][0] - vals[0][0]
        ambig_counts[diff_val] += 1
        if diff_val > ambig_threshold:
            species = vals[0][1]
            best_read = dbs.to_be_resolved_reads[read.qname][species]
            # assigned.write(best_read)
            specific_files[species].write(best_read)
            species_counts[species] += 1
        else:
            species_counts['ambig'] += 1
            best = vals[0][0]
            ambig_names = []
            for val, spec in vals:
                if val - best > ambig_threshold:
                    break
                if spec not in ambig_names:
                    ambig_names.append(spec)

            ambig_names = tuple(ambig_names)
            ambig_types[ambig_names] += 1
            for amb_read in dbs.to_be_resolved_reads[read.qname].itervalues():
                ambig.write(amb_read)

    # Clean up the dictionaries
    dbs.to_be_resolved_vals.pop(read.qname)
    dbs.to_be_resolved_counts.pop(read.qname)
    dbs.to_be_resolved_reads.pop(read.qname)

ambig_threshold = 3

for fname in sys.argv[1:]:
    print fname
    samfile = pysam.Samfile(fname, 'rb')
    references = samfile.references
    dir = path.dirname(fname)
    # assigned = pysam.Samfile(path.join(dir, 'assigned.bam'), 'wb',
    # template=samfile)
    ambig = pysam.Samfile(path.join(dir, 'ambiguous.bam'), 'wb',
                          template=samfile)
    specific_files = my_defaultdict(pysam.Samfile,
                                    path.join(dir, 'assigned_%s.bam'),
                                    {'template': samfile,
                                     'mode': 'wb'})

    to_be_resolved_reads = defaultdict(dict)
    to_be_resolved_vals = defaultdict(lambda: defaultdict(lambda: 1000))
    to_be_resolved_counts = Counter()
    to_be_resolved_reads2 = defaultdict(dict)
    to_be_resolved_vals2 = defaultdict(lambda: defaultdict(lambda: 1000))
    to_be_resolved_counts2 = Counter()
    species_counts = Counter()
    ambig_counts = Counter()
    ambig_types = Counter()

    print "Measuring file size"
    start = samfile.tell()
    maxval = path.getsize(fname) * 2**16  # I don't know why it's off by 2^16
    pbar = ProgressBar(maxval=maxval - start + 2**16,
                       widgets=[fname, ': ', Percentage(), ' ', Bar(), ' ',
                                ETA(), ' '])
    pbar.start()

    for read in samfile:
        pbar.update(samfile.tell() - start)
        process_read(read)
    pbar.finish()
    print
    print "Species assignments in %s: %s" % (fname, species_counts)
    print "Ambiguity distribution: ", ambig_counts
    print "Ambiguity types: ", ambig_types.most_common(50)

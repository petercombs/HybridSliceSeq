from __future__ import print_function
from collections import defaultdict
from pysam import Samfile
from progressbar import ProgressBar as PBar
import pandas as pd
import numpy as np
from argparse import ArgumentParser
from sys import stdout

def parse_annotation(annot):
    retval = {}
    for i in annot.strip(';').split(';'):
        key, val = i.strip().split(' ', maxsplit=1)
        retval[key.strip()] = val.strip().strip(';"')
    return retval


def get_supported_exons(read, exons_in_gene):
    exon_num = 0
    exons_iter = iter(exons_in_gene)
    supported_exons = []
    supported_lens = []
    excluded_exons = []
    cur_pos = read.pos
    cur_exon = next(exons_iter)
    try:
        while cur_exon[1] < cur_pos:
            #print("Skipping", cur_exon, cur_pos)
            cur_exon = next(exons_iter)
            #print("Skipped, now on:", cur_exon, cur_exon[1] < cur_pos)
            exon_num += 1
        #print(read.cigartuples)
        for cigartype, length in read.cigartuples:
            # Match = 0, Deletion = 2, Skipped = 3
            #print(cigartype, length, cur_exon, cur_pos,
                  #cur_pos <= cur_exon[1] <= cur_pos+length,
                  #cur_pos < cur_exon[1] < cur_pos + length,
                 #)
            if cigartype == 0:
                while ((cur_pos <= cur_exon[1] <= cur_pos+length) or
                       (cur_pos <= cur_exon[0] <= cur_pos+length) or
                       (cur_exon[0] < cur_pos < cur_exon[1])):
                    supported_lens.append(cur_exon[1] - cur_exon[0])
                    supported_exons.append(cur_exon[2])
                    #print("Adding", cur_exon)
                    cur_exon = next(exons_iter)
                cur_pos += length
            if cigartype==3 or cigartype == 2:
                while cur_pos < cur_exon[1] < cur_pos + length:
                    #print("Excluding", cur_exon)
                    excluded_exons.append(cur_exon[2])
                    cur_exon = next(exons_iter)
                cur_pos += length
    except StopIteration:
        #print("StoppedIteration", cur_exon)
        pass
    finally:
        return (
            supported_lens,
            supported_exons,
            (excluded_exons if supported_exons else [])
        )


def get_exon_dictionary(gtf_file):
    current_gene = ''
    all_exons_by_chrom_by_gene = defaultdict(list)
    exon_names = {}

    num_exons = 0
    for line in open(gtf_file):
        data = line.strip().split('\t')
        left = int(data[3])
        right = int(data[4])
        annot = parse_annotation(data[-1])
        if data[2] == 'aggregate_gene':
            current_exons = []
            current_gene = annot['gene_id']
            # Note here that what *should* happen is due to python's
            # pass-by-reference magic, we should be able to keep editing
            # current_exons even after it's been inserted into the dictionary
            all_exons_by_chrom_by_gene[data[0]].append(
                (left, right, current_gene, current_exons)
            )
        elif data[2] == 'exonic_part':
            current_exons.append((left, right, num_exons,
                                  '{}_{}'.format(
                                      annot['gene_id'], annot['exonic_part_number']
                                  )
                                 ))
            exon_names[num_exons] = '{}_{}'.format( annot['gene_id'], annot['exonic_part_number'])
            num_exons += 1
    return all_exons_by_chrom_by_gene, exon_names


def parse_overlapping_reads(samfile, raw_exon_counts, exons_on_chrom):
    for read in samfile:
        chrom = read.reference_name
        for read in samfile:
            for gene_left, gene_right, gene_name, exons in exons_on_chrom[chrom]:
                positions = read.positions
                if gene_left <= positions[0] < positions[-1] <= gene_right:
                    old_exon_right = -1
                    for exon_left, exon_right in exons:
                        if exon_right < positions[0] or positions[-1] < exon_left:
                            continue
                        elif exon_left in positions or exon_right in positions:
                            raw_exon_counts.ix[(chrom, exon_left, exon_right,
                                                gene), 'INCLUDED'] += 1
                        else:
                            pass

def iterate_over_samfile(samfile, exon_dictionary, exon_names, timeit=1e100):
    samfile.reset()
    curr_reference = ''
    all_exons = pd.DataFrame(index=exon_names.keys(),
                             columns=['SUPPORTED', 'EXCLUDED',
                                     'N_SUPPORTED', 'N_EXCLUDED'],
                             data=0.0, dtype=float,
                             #data={'SUPPORTED': 0.0, 'EXCLUDED': 0.0,
                                   #'N_SUPPORTED': 0, 'N_EXCLUDED': 0}
                            )
    num_exons = len(exon_names)
    arr_supported = np.zeros(num_exons, dtype=float)
    arr_excluded = np.zeros(num_exons, dtype=float)
    arr_n_supported = np.zeros(num_exons, dtype=int)
    arr_n_excluded = np.zeros(num_exons, dtype=int)
    references=samfile.references
    pb = PBar(maxval=samfile.mapped)
    pb.start()
    #pb = lambda x: x
    for i, read in enumerate(samfile):
        if i % 1e5 == 0:
            pb.update(i)
        refid = read.reference_id
        if refid != curr_reference:
            curr_reference = refid
            genes_on_chrom = exon_dictionary[references[refid]]
            current_position_in_annotation = 0
        curpos = read.pos
        while current_position_in_annotation < len(genes_on_chrom) and genes_on_chrom[current_position_in_annotation][1] < curpos:
            current_position_in_annotation += 1
        for left, right, genedata, exons in genes_on_chrom[current_position_in_annotation:]:
            if curpos < left:
                break
            readlen = sum(i for t, i in read.cigartuples if t == 0)
            lens, supported, excluded = get_supported_exons(read, exons)
            for length, ix in zip(lens, supported):
                arr_supported[ix] += 1/(length + readlen - 1)
                arr_n_supported[ix] += 1
            for ix in excluded:
                arr_excluded[ix] += 1/(readlen - 1)
                arr_n_excluded[ix] += 1
        if i > timeit:
            break
    pb.finish()
    all_exons['SUPPORTED'] = arr_supported
    all_exons['EXCLUDED'] = arr_excluded
    all_exons['N_SUPPORTED'] = arr_n_supported
    all_exons['N_EXCLUDED'] = arr_n_excluded
    all_exons.rename(index=exon_names, inplace=True)
    return all_exons

def iterate_over_references(samfile, exon_dictionary):
    samfile.reset()
    references = samfile.references
    all_exons = pd.DataFrame(index=[exon
                                    for genes in exon_dictionary.values()
                                    for _, _, _, exons in genes
                                    for _, _, exon in exons
                                   ],
                             columns=['SUPPORTED', 'EXCLUDED', 'N_SUPPORTED',
                                      'N_EXCLUDED'],
                             data=0)

    for reference in PBar()(references):
        genes_on_chrom = exon_dictionary[reference]
        current_position_in_annotation = 0
        for read in samfile.fetch(reference):
            curpos = read.pos
            while current_position_in_annotation < len(genes_on_chrom) and genes_on_chrom[current_position_in_annotation][1] < curpos:
                current_position_in_annotation += 1
            for left, right, genedata, exons in genes_on_chrom[current_position_in_annotation:]:
                if curpos < left:
                    break
                supported, excluded = get_supported_exons(read, exons)
                all_exons.ix[supported, 'SUPPORTED'] += 1
                all_exons.ix[excluded, 'EXCLUDED'] += 1


    return

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('samfile', type=Samfile)
    parser.add_argument('gtf_file')
    parser.add_argument('--outfile', '-o', default=stdout)
    parser.add_argument('--min-reads', default=20, type=int)
    parser.add_argument('--timeit', default=1e100, type=float)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    exon_dict, exon_names = get_exon_dictionary(args.gtf_file)
    all_exons = iterate_over_samfile(args.samfile, exon_dict, exon_names, args.timeit)
    all_exons.index.name = 'exon_id'
    all_exons['psi'] = all_exons.SUPPORTED / (all_exons.SUPPORTED +
                                              all_exons.EXCLUDED)
    all_exons.ix[(all_exons.N_SUPPORTED + all_exons.N_EXCLUDED) < args.min_reads,
                 'psi'] = pd.np.nan
    all_exons.to_csv(args.outfile, sep='\t', na_rep='NA')






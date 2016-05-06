from pysam import Samfile
from sys import argv
from collections import defaultdict
from os import path
import svgwrite as svg

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
    if snpfile.endswith('.bed'):
        for line in open(snpfile):
            chrom, _, start, refalt = line.strip().split()
            if true_hets.get((chrom, start), True):
                snps[chrom][int(start)] = refalt.split('|')

    return snps

read_height = 5
x_scale = 1
spacing = 5
pos_color = 'blue'
neg_color = 'red'
unk_color = 'gray'
max_rows = 100


def insert_reads(read1, read2, reads_by_phase):
    if read2 is not None and read2.pos < read1.pos:
        read1, read2 = read2, read1
    for i, row in enumerate(reads_by_phase):
        if read1.reference_start > row[-1].reference_end + spacing:
            if read2 is None:
                row.append(read1)
            elif read1.pos < read2.pos:
                row.append(read1)
                row.append(read2)
            else:
                row.append(read2)
                row.append(read1)

            break
        elif i > max_rows:
            row.append(read1)
            if read2 is not None:
                row.append(read2)
            break
    else:
        reads_by_phase.append([read1])

def draw_read(read, dwg, x_start_coord, y_coord, phase_color, snps=None,
              with_snps_only=False, last_read=None):
    blocks = read.blocks
    slice_start, slice_end = 0, len(blocks)
    snp_locs = []
    for read_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
        if snps is not None and (ref_pos+1) in snps:
            snp_locs.append((read_pos, ref_pos))
    if with_snps_only and len(snp_locs) == 0:
        return
    if read.is_reverse:
        slice_start += 1
        block_start, block_end = blocks[0]
        last_stop = block_end
        dwg.add(dwg.polygon(
            [(x_scale * (block_start - x_start_coord), y_coord + read_height/2),
             (x_scale * (min(block_start+read_height/2, block_end) - x_start_coord), y_coord+read_height),
             (x_scale * (block_end - x_start_coord), y_coord+read_height),
             (x_scale * (block_end - x_start_coord), y_coord),
             (x_scale * (min(block_start+read_height/2, block_end) - x_start_coord), y_coord),
             (x_scale * (block_start - x_start_coord), y_coord+read_height/2)
            ],
            style='stroke-width=0;fill:'+phase_color,
            id=read.qname,
        ))
    else:
        slice_end -= 1
        block_start, block_end = blocks[-1]
        dwg.add(dwg.polygon(
            [
                (x_scale * (block_start - x_start_coord), y_coord + read_height),
                (x_scale * (max(block_start, block_end - read_height/2) - x_start_coord), y_coord+read_height),
                (x_scale * (block_end - x_start_coord), y_coord+read_height/2),
                (x_scale * (max(block_start, block_end - read_height/2) - x_start_coord), y_coord),
                (x_scale * (block_start - x_start_coord), y_coord),
                (x_scale * (block_start - x_start_coord), y_coord+read_height),
            ],
            style='stroke-width=0;fill:'+phase_color,
        ))
        last_stop = blocks[0][0]

    for i in range(slice_start, slice_end):
        block_start, block_stop = blocks[i]
        dwg.add(dwg.line((x_scale*(last_stop - x_start_coord),
                          y_coord + 0.5 * read_height),
                         (x_scale*(block_start-x_start_coord),
                          y_coord + 0.5 * read_height),
                         style='stroke-width=2;stroke:'+phase_color,
                         id=read.qname,
                        ))
        dwg.add(dwg.rect((x_scale*(block_start - x_start_coord),
                          y_coord),
                         (x_scale*(block_stop - block_start),
                          read_height),
                         style='stroke-width=0;fill:'+phase_color,
                         id=read.qname,
                        ))
    if snps:
        for read_pos, ref_pos in snp_locs:
            alleles = snps[ref_pos + 1]
            if read.seq[read_pos] in alleles:
                snp_color = [neg_color, pos_color][alleles.index(read.seq[read_pos])]
            else:
                snp_color = 'black'
            dwg.add(dwg.line(
                (x_scale * (ref_pos - x_start_coord), y_coord),
                (x_scale * (ref_pos - x_start_coord), y_coord + read_height),
                style='stroke-width:1; stroke:' + snp_color,
            ))

    if last_read is not None and read.qname == last_read.qname:
        dwg.add(dwg.line((x_scale * (last_read.reference_end - x_start_coord),
                          y_coord + read_height/2),
                         (x_scale * (read.reference_start - x_start_coord),
                          y_coord + read_height/2),
                         style='stroke-width:1;stroke:black',
                        ))



if __name__ == "__main__":
    phase_pos = []
    phase_neg = []
    phase_unk  = []

    # Note that phases are 0, 1, and -1
    phase_all = [phase_unk, phase_pos, phase_neg]
    snps = get_snps(argv[1])

    start_coord = 1e99

    unmatched_reads = [{}, {}]
    # Note that in order to keep track of the un-phased reads, we need to do a
    # 2-pass approach to know the height of all of the classes
    gene_chrom = 'dmel_2L'
    gene_coords = (8258373,8301079)
    num_reads = 0
    for read in (Samfile(argv[2])
                 .fetch(gene_chrom, gene_coords[0], gene_coords[1])):
        if (not ((gene_coords[0] <= read.reference_start <= gene_coords[1])
                 and (gene_coords[0] <= read.reference_end <= gene_coords[1]))):
            continue
        phase = get_phase(read, snps[gene_chrom])

        # read.is_read1 = 1 if read_1, 0 if read_2
        if read.qname in unmatched_reads[read.is_read1]:
            other_read, other_phase = unmatched_reads[read.is_read1].pop(read.qname)

            if phase is None and other_phase is None:
                continue
            elif other_phase is None:
                pass
            elif phase is None:
                phase = other_phase
            elif phase != other_phase:
                phase = 0
            else:
                pass

            reads_by_phase = phase_all[phase]
            insert_reads(read, other_read, reads_by_phase)
        else:
            unmatched_reads[read.is_read2][read.qname] = (
                read, phase
            )


        if phase is None:
            continue

        start_coord = min(start_coord, read.reference_start)

    for reads in unmatched_reads:
        for qname in reads:
            read, phase = reads[qname]
            if phase is None:
                continue
            reads_by_phase = phase_all[phase]
            insert_reads(read, None, reads_by_phase)




    max_depth_pos = len(phase_pos)
    max_depth_neg = len(phase_neg)
    max_depth_unk = len(phase_unk)

    dwg = svg.Drawing(argv[2].replace('.bam', '_phased.svg'))

    y_start = 10 + max_depth_neg * 1.2*read_height

    last_read = None
    for row_num, row in enumerate(phase_neg):
        for read in row:
            if read is None:
                continue
            draw_read(read, dwg,
                      start_coord, y_start - row_num * 1.2 * read_height,
                      neg_color,
                      last_read = last_read
                     )
            last_read = read

    y_start += 3 * read_height
    for row_num, row in enumerate(phase_unk):
        for read in row:
            draw_read(read, dwg,
                      start_coord, y_start + row_num * 1.2 * read_height,
                      unk_color,
                      snps=snps[gene_chrom],
                      last_read = last_read
                     )
            last_read = read
    y_start += max_depth_unk * 1.2 * read_height + 3 * read_height

    for row_num, row in enumerate(phase_pos):
        for read in row:
            draw_read(read, dwg,
                      start_coord, y_start + row_num * 1.2 * read_height,
                      pos_color,
                      last_read = last_read
                     )
            last_read = read

    dwg.save()



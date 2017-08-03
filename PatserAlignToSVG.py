from __future__ import print_function
import svgwrite as svg
from Bio import  AlignIO, SeqIO
from argparse import ArgumentParser
from glob import glob
from os import path
from numpy import where, argwhere, diff, zeros, random, ceil
import pandas as pd
from progressbar import ProgressBar
import itertools as it
from collections import defaultdict, Counter
from pprint import pprint

seq_colors = defaultdict(
    lambda : 'gray',
    {'T' : 'red', 'A': 'green', 'C': 'blue', 'G':'gold'}
)

def get_header_size(filename):
    for i, line in enumerate(open(filename)):
        if 'position=' in line:
            return i
    return 0

def has_match(pos, patser1, patser2, pos1, pos2, dist=5, reltol=.2):
    score = patser1.score[pos]
    try:
        new_pos = pos2[argwhere(pos1 == pos)][0,0]
    except IndexError:
        return False
    return bool(len(patser2
            .query('({lp} < pos) and (pos < {hp}) and ({ls} < score) and (score < {hs})'
                .format(
                    lp=new_pos - dist,
                    hp=new_pos + dist,
                    ls=(1-reltol)*score,
                    hs=(1+reltol)*score,
                    )
                )))



def parse_best_aligns(aligns_by_seq1, take_best=min):
    out = []
    for seq1 in aligns_by_seq1:
        #print(seq1, [a[:-1] + (a[-1][1].id, ) for a in aligns_by_seq1[seq1]])
        bestal = take_best(aligns_by_seq1[seq1])
        bestal = bestal[-1]
        bestlen = len(bestal[0])
        pos1 = zeros(bestlen, dtype=int)
        pos2 = zeros(bestlen, dtype=int)
        i1 = 0
        i2 = 0
        for i, (b1, b2) in enumerate(zip(bestal[0], bestal[1])):
            pos1[i] = i1
            pos2[i] = i2
            i1 += b1 != '-'
            i2 += b2 != '-'

        out.append(((bestal[0].id, pos1), (bestal[1].id, pos2)))
    return out





def parse_args():
    store_true = 'store_true'
    parser = ArgumentParser()
    data_opts = parser.add_argument_group('Data file options')
    plot_opts = parser.add_argument_group('Plot options')
    parser.add_argument('--tf', '-t', default=[], nargs='*')
    parser.add_argument('--tf-names', '-T', default=[], nargs='*')
    parser.add_argument('--comp1', '-1', default='anterior')
    parser.add_argument('--comp2', '-2', default='posterior')
    data_opts.add_argument('--sequence', '-s', default=False, action=store_true,
                           help='Patser output the sequence of the match as well')
    data_opts.add_argument('--fasta', '-F', default=None)
    data_opts.add_argument('--needleall', '-N', default=False, action=store_true,
                           help='Input alignments are actually needleall output'
                           " files in EMBOSS srs format")
    data_opts.add_argument('--meme-suite', '-M', default=False,
                           action=store_true,
                           help='Directory is actually the output of FIMO or '
                           'another MEME suite program')
    data_opts.add_argument('--bed-track', '-b', default=False)
    data_opts.add_argument('--gtf-track', '-g', default=False)
    data_opts.add_argument('--coordinates-bed', default=False)
    plot_opts.add_argument('--x-scale', '-x', type=float, default=.1)
    plot_opts.add_argument('--y-scale', '-y', type=float, default=1.0)
    plot_opts.add_argument('--y-sep', '-Y', type=float, default=50)
    plot_opts.add_argument('--x-ticks', '-X', type=float, default=1000,
                           help='Put a marker every X-TICKS bases on the'
                           'upper-most strand')
    plot_opts.add_argument('--margin', '-m', type=float, default=10)
    plot_opts.add_argument('--bar-width', '-w', type=float, default=1)
    plot_opts.add_argument('--draw-alignment', '-A', default=False, action=store_true)
    plot_opts.add_argument('--no-draw-binding', dest='draw_binding', default=True, action='store_false')
    plot_opts.add_argument('--match-dim', '-d', type=float, default=0.5)
    plot_opts.add_argument('--rescale-bars', '-R', default=False,
                           action=store_true)
    plot_opts.add_argument('--background', '-B', default=False,
                           action=store_true)
    plot_opts.add_argument('--bed-height', '-H', default=0, type=float)
    plot_opts.add_argument('--make-png', '-P', default=False, action=store_true)
    parser.add_argument('--print-argv', '-v', default=False,
                        action=store_true,)
    parser.add_argument('--num-tf-cols', '-C', default=3, type=int,
                        help="Number of columns in TF labels")
    parser.add_argument('alignments', type=open,
                        help='File containing alignments in MUMmer .delta format')
    parser.add_argument('patser_directory', type=str)


    args =  parser.parse_args()
    if args.print_argv:
        from sys import argv
        print(" ".join(argv))

    if args.tf == []:
        args.tf = ['bcd', 'cad', 'gt', 'hb', 'kni', 'kr', 'hkb', 'tll', 'D', 'ftz', 'h', 'prd', 'run', 'slp1']

    if args.fasta:
        args.fasta = {rec.id: rec for rec in SeqIO.parse(args.fasta, 'fasta')}

    if args.tf_names == []:
        args.tf_names = args.tf
    elif len(args.tf_names) != len(args.tf):
        raise ValueError("Different number of TF labels and TFs: {} and {}"
                         .format(args.tf_names, args.tfs))

    if args.coordinates_bed:
        args.has_coordinates = True
        args.coordinates_bed = pd.read_table(
            args.coordinates_bed,
            header=None, names=['chr', 'start', 'stop', 'feature', 'score',
                                'strand', 'thickStart', 'thickEnd'],
            index_col='feature',
        )
    else:
        args.has_coordinates = False


    if args.has_coordinates and args.bed_track:
        args.bed_track = pd.read_table(
            args.bed_track,
            comment='#',
            header=None, names=['chr', 'start', 'stop', 'name', 'score', ],
        )
        args.draw_bed = True
    else:
        args.draw_bed = False

    if args.has_coordinates and args.gtf_track:
        args.gtf_track = pd.read_table(
            args.gtf_track,
            header=None, names=['chr', 'source', 'type', 'start', 'stop',
                                'score', 'strand', 'offset', 'annotation']
        )
        args.draw_gtf = True
    else:
        args.draw_gtf = False


    return args

def get_drawing_size(parsed_arguments, alignments):

    width = 200
    height = 0

    width += 2* parsed_arguments.margin
    width += args.x_scale * max(len(a[0][1]) for a in alignments)
    width = max(width, 500)

    height += .5
    height += len(alignments)
    height += len(alignments) * .5 * parsed_arguments.draw_alignment
    height += (0.5 * ceil(len(parsed_arguments.tf)/3))
    height += len(alignments) * .9 * args.draw_bed
    height += len(alignments) * .9 * args.draw_gtf

    height *= 1.5 * parsed_arguments.y_sep
    height += 200

    return width, height

def draw_tf_cols(dwg, dwg_groups, args, y_start):
    n_cols = args.num_tf_cols
    x_start = args.margin
    delta_y = 1.5 * args.y_sep

    for i, (tf, tf_name) in enumerate(zip(args.tf, args.tf_names)):
        dwg_groups[tf].add(dwg.rect(
            (x_start + 120  * (i % n_cols), y_start),
            (.3*delta_y, .3*delta_y),
            **{'fill' : colors[i],
                #'stroke-width': "0",
                }
            #style='fill:{}; stroke-width:0;'.format(colors[i])

            ))
        dwg_groups[tf].add(dwg.text(tf_name,
            (x_start + 10 + 120  * (i % n_cols)+.3*delta_y, y_start+0.3*delta_y)
            ))
        y_start += delta_y * ((i%n_cols) == (n_cols-1))

    for tf in dwg_groups:
        dwg.add(dwg_groups[tf])
    for tf in dwg_groups:
        grp = dwg_groups[tf].copy()
        grp.attribs['class'] += ' hover_group'
        dwg.add(grp)
    for line in lines:
        dwg.add(line)

    return y_start


def draw_gtf_track(args, dwg, n1, n2, pos1, pos2, y_start):
    g = dwg.g()
    minus_strand = False
    '''
    print(args.coordinates_bed.index,
          n1, n2,
          args.coordinates_bed.ix[:, 'strand'],
          file=stderr)
          '''
    if n1 in args.coordinates_bed.index:
        if args.coordinates_bed.ix[n1, 'strand'] == '-':
            bed_pos = args.coordinates_bed.ix[n1, 'stop'] - pos1
            minus_strand = True
        else:
            bed_pos = pos1 + args.coordinates_bed.ix[n1, 'start']
    elif n2 in args.coordinates_bed.index:
        if args.coordinates_bed.ix[n2, 'strand'] == '-':
            bed_pos = args.coordinates_bed.ix[n2, 'stop'] - pos2
            minus_strand = True
        else:
            bed_pos = pos2 + args.coordinates_bed.ix[n2, 'start']
    else:
        raise RuntimeError("Cannot find coordinates for: {}".format(n1))


    #print(bed_pos, file=stderr)
    #print(args.coordinates_bed.to_csv(sep='\t'))
    '''
    print(pos1 + args.coordinates_bed.ix[n2, 'start'], file=stderr)
    print(pos2 + args.coordinates_bed.ix[n2, 'start'], file=stderr)
    '''
    height = 0.15 * delta_y
    for ix, row in args.gtf_track.iterrows():
        #print(minus_strand, row.strand, row.annotation, file=stderr)
        if minus_strand:
            left = min(row.stop, bed_pos[0])
            right = max(row.start, bed_pos[-1])
            x_starts = where(bed_pos == left)[0]
            x_ends = where(bed_pos == right)[0]
            x_start = x_starts[0] if len(x_starts) else 0
            x_end = x_ends[-1] if len(x_ends) else len(bed_pos)
            width = x_end - x_start
            if row.strand == '+':
                x_start = x_start * args.x_scale + args.margin
                x_end = x_start + width * args.x_scale
                r = dwg.path(
                    'M{r},{t} H{rm} L{l},{m} L{rm},{b} H{r} Z'
                    .format(l=x_start, r = x_end, rm=x_start + height/2,
                            t=y_start-height, m=y_start - height/2,
                            b=y_start),
                    **{'fill': 'gray', 'class_': 'hoverstroke'}
                )
            else:
                x_start = x_start * args.x_scale + args.margin
                x_end = x_start + width * args.x_scale
                r = dwg.path(
                    'M{l},{t} H{rm} L{r},{m} L{rm},{b} H{l} Z'
                    .format(l=x_start, r = x_end, rm=x_end-height/2,
                            t=y_start-height, m=y_start - height/2,
                            b=y_start),
                    **{'fill': 'gray', 'class_': 'hoverstroke',
                      }
                )
        else:
            left = max(row.start, bed_pos[0])
            right = min(row.stop, bed_pos[-1])
            x_starts = where(bed_pos == left)[0]
            x_ends = where(bed_pos == right)[0]
            x_start = x_starts[0] if len(x_starts) else len(bed_pos)
            x_end = x_ends[-1] if len(x_ends) else 0
            #print(left, right, '---', x_start, x_end, file=stderr)
            width = x_end - x_start
            if row.strand == '+':
                x_start = x_start * args.x_scale + args.margin
                x_end = x_start + width * args.x_scale
                r = dwg.path(
                    'M{l},{t} H{rm} L{r},{m} L{rm},{b} H{l} Z'
                    .format(l=x_start, r = x_end, rm=x_end - height/2,
                            t=y_start-height, m=y_start - height/2,
                            b=y_start),
                    **{'fill': 'gray', 'class_': 'hoverstroke'}
                )
            else:
                x_start = x_start * args.x_scale + args.margin
                x_end = x_start + width * args.x_scale
                r = dwg.path(
                    'M{r},{t} H{rm} L{l},{m} L{rm},{b} H{r} Z'
                    .format(l=x_start, r = x_end, rm=x_start + height/2,
                            t=y_start-height, m=y_start - height/2,
                            b=y_start),
                    **{'fill': 'gray', 'class_': 'hoverstroke',
                      }
                )
        r.add(svg.base.Title(row.annotation))
        g.add(r)

    y_start += args.y_sep * .9
    dwg.add(g)


    return y_start


def draw_bed_track(args, dwg, n1, n2, pos1, pos2, y_start):
    minus_strand = False
    if n1 in args.coordinates_bed.index:
        chrom = args.coordinates_bed.ix[n1, 'chr']
        if args.coordinates_bed.ix[n1, 'strand'] == '-':
            bed_pos = args.coordinates_bed.ix[n1, 'stop'] - pos1
            minus_strand = True
        else:
            bed_pos = pos1 + args.coordinates_bed.ix[n1, 'start']
    elif n2 in args.coordinates_bed.index:
        chrom = args.coordinates_bed.ix[n2, 'chr']
        if args.coordinates_bed.ix[n2, 'strand'] == '-':
            bed_pos = args.coordinates_bed.ix[n2, 'stop'] - pos1
            minus_strand = True
        else:
            bed_pos = pos1 + args.coordinates_bed.ix[n2, 'start']
    else:
        raise RuntimeError("Cannot find coordinates for: {}".format(n1))

    bed_track_lo = args.bed_track.query(
        ('chr == "{chrom}" and start <= {pos} and {pos} < stop')
        .format(chrom=chrom, pos=min(bed_pos))
    )

    bed_track_hi = args.bed_track.query(
        ('chr == "{chrom}" and start < {pos} and {pos} <= stop')
        .format(chrom=chrom, pos=max(bed_pos))
    )

    bed_track = args.bed_track.ix[bed_track_lo.index[0]:bed_track_hi.index[-1]]
    if args.bed_height:
        bed_track_score_norm = args.bed_height
    else:
        bed_track_score_norm = max(bed_track.score)

    for ix, row in bed_track.iterrows():
        if minus_strand:
            pass
        else:
            x_starts = where(bed_pos == row.start)[0]
            x_ends = where(bed_pos == row.stop)[0]
            x_start = x_starts[0] if len(x_starts) else 0
            x_end = x_ends[-1] if len(x_ends) else len(bed_pos)

            width = (x_end - x_start)
            height = row.score / bed_track_score_norm * .4 * args.y_sep
            dwg.add(dwg.rect(
                (x_start*args.x_scale + args.margin, y_start-height),
                (width*args.x_scale, height),
                **{'fill': 'gray'}
            ))

    y_start += args.y_sep * .9


    return y_start

def draw_alignments(args, dwg, n1, n2, pos1, pos2, y_start):
    x_start = args.margin
    fasta1 = path.join(args.patser_directory, n1 + '.fasta')
    fasta2 = path.join(args.patser_directory, n2+'.fasta')
    if path.exists(fasta1) and path.exists(fasta2):
        seq1 = SeqIO.read(fasta1, 'fasta').seq
        seq2 = SeqIO.read(fasta2, 'fasta').seq
    elif args.fasta:
        if n1 in args.fasta and n2 in args.fasta:
            seq1 = args.fasta[n1].seq
            seq2 = args.fasta[n2].seq
        else:
            print("Can't find both sequences: {} and {}".format(n1, n2))
            seq1 = defaultdict(lambda : 'N')
            seq2 = defaultdict(lambda : 'N')
    else:
        seq1 = defaultdict(lambda : 'N')
        seq2 = defaultdict(lambda : 'N')

    lines.append(dwg.line(
        (x_start, y_start),
        (x_start + args.x_scale * len(pos1), y_start)))

    indels = diff(pos2) - diff(pos1)
    id_start = 0
    # Draw Indels
    for val, group in it.groupby(indels):
        for i, _ in enumerate(group):
            pass
        i += 1
        dwg.add(dwg.rect(
            (x_start + x_scale * id_start, y_start - .1 * args.y_sep * (val < 0)),
            (x_scale * i, .1 * args.y_sep * (val != 0)),
            fill="grey"
            ))
        if val == 0:
            for j in range(id_start, id_start + i):
                if str(seq1[pos1[j]]) != str(seq2[pos2[j]]):
                    dwg.add(dwg.line(
                        (x_start + x_scale * j, y_start - .3 * args.y_sep),
                        (x_start + x_scale * j, y_start),
                        style="stroke-width:1; stroke:{};".format(seq_colors[str(seq1[pos1[j]])]),
                        ))
                    dwg.add(dwg.line(
                        (x_start + x_scale * j, y_start),
                        (x_start + x_scale * j, y_start + .3 * args.y_sep),
                        id='{}:{}--{}>{}'.format(pos1[j], pos2[j], seq1[pos1[j]], seq2[pos2[j]],),
                        style="stroke-width:1; stroke:{};".format(seq_colors[str(seq2[pos2[j]])]),
                        ))

        id_start += i
    y_start += 0.5 * delta_y
    return y_start

if __name__ == "__main__":
    args = parse_args()

    if args.meme_suite:
        meme1 = (
            pd.read_table( path.join(args.patser_directory, args.comp1, 'fimo.txt'),)
            .rename(columns=lambda x: x.strip('#'))
            .rename(columns=lambda x: x.replace(' ', '_'))
        )
        meme1['pos'] = (meme1.start + meme1.stop)//2
        meme1['tf'] = meme1['pattern_name']
        meme2 = (
            pd.read_table( path.join(args.patser_directory, args.comp2, 'fimo.txt'),)
            .rename(columns=lambda x: x.strip('#'))
            .rename(columns=lambda x: x.replace(' ', '_'))
        )
        meme2['pos'] = (meme2.start + meme2.stop)//2
        meme2['tf'] = meme2['pattern_name']



    #header_size = get_header_size(glob(path.join(args.patser_directory, '*'))[0])
    patser_args = dict(
        index_col=['pos'],
        names=['seq', '_', 'pos', '__', 'score', '___', 'pval']
        + (['____', 'bases'] if args.sequence else [])
        ,
                usecols=[0,2,4,6]+([8] if args.sequence else []),
                #skiprows=header_size,
                na_values=[' ', '', '-', '='],
                sep=' +',
                engine='python',
                #dtype={'pos':int, 'score':float, 'pval':float},
            )


    x_scale = args.x_scale
    y_scale = args.y_scale
    x_start = args.margin
    y_start = 1.5 * args.y_sep
    delta_y = 1.5 * args.y_sep
    match_dim = args.match_dim

    colors=[
            'darkorange',
            'blue',
            'gold',
            'aqua',
            'firebrick',
            'darkblue',
            'forestgreen',
            'lime',
            'tomato',
            'chartreuse',
            'darkslateblue',
            'brown',
            'crimson',
            'burlywood',
            'mediumspringgreen',
            'mediumorchid',
            'darkcyan',
            'black',
            ]

    alignments = []
    if args.needleall:
        aligns_by_seq1 = defaultdict(list)
        aligns_by_seq2 = defaultdict(list)
        all_aligns = {}
        for align in AlignIO.parse(args.alignments, 'emboss'):
            all_aligns[(align[0].id, align[1].id)] = align

        a1_name = ''
        a2_name = ''
        args.alignments.seek(0)
        for line in args.alignments:
            if line.startswith('# 1:'):
                a1_name = line.strip().split(':')[1].strip()
            elif line.startswith('# 2:'):
                a2_name = line.strip().split(':')[1].strip()
            elif line.startswith('# Score:'):
                align = all_aligns[a1_name, a2_name]
                aligns_by_seq1[a1_name].append((
                    -(Counter(align[0])['-']+Counter(align[1])['-']),
                    -len(align[0]),
                    float(line.strip().split(':')[1]),
                    random.rand(),
                    align,
                ))
                aligns_by_seq2[a2_name].append((
                    sum(align[0][i] == align[1][i]
                        for i in range(len(align[0])))
                    -sum((align[0][i] == '-') or (align[1][i] == '-')
                        for i in range(len(align[0]))),
                    -(Counter(align[0])['-']+Counter(align[1])['-']),
                    -len(align[0]),
                    float(line.strip().split(':')[1]),
                    random.rand(),
                    align,
                ))
        alignments = parse_best_aligns(aligns_by_seq2, max)
        n_aligns = len(alignments)
    else:
        try:
            from ParseDelta import parse_records
        except err as ImportError:
            print("Ask Peter Combs for the ParseDelta.py file")
            raise err

        for n_aligns, a in enumerate(parse_records(args.alignments)):
            alignments.append(a)
    pb = ProgressBar(maxval=(n_aligns+1)*len(args.tf))
    print("Number of TFs", len(args.tf))

    width, height = get_drawing_size(args, alignments)
    dwg = svg.Drawing(
        args.alignments.name+'.svg',
        size=(width, height),
    )
    dwg_groups = {}
    for tf in args.tf:
        dwg_groups[tf] = dwg.g(class_=tf)
    dwg.add(dwg.style(
        'line{stroke-width:1;stroke:#000000;}\n'
        '.hover_group{opacity:0.5;} \n'
        '.hover_group:hover \n'
        '{\n\topacity:1;\n'
        '\tstroke-width:1!important;'
        '\n\tstroke:#000001;\n}\n'
        '.hoverstroke{stroke-width:0; stroke:black; opacity:.5} \n'
        '.hoverstroke:hover{stroke-width:1; opacity:1}\n'
        '.fixed{stroke:green;}\n'
        '.segregating{stroke:red;opacity:0.01}\n'
        '.diverged{stroke:black;opacity:0.1}\n'
        + ('.matched{{opacity: {} }}\n'.format(args.match_dim))
    ))

    prog = 0
    lines = []
    tf_changes = {}
    for ((n1, pos1), (n2, pos2)) in sorted(alignments,
                                           key=lambda x: (x[0][0], x[1][0])):
        if args.draw_alignment:
            y_start = draw_alignments(args, dwg, n1, n2, pos1, pos2, y_start)

        if args.draw_bed:
            y_start = draw_bed_track(args, dwg, n1, n2, pos1, pos2, y_start)

        if args.draw_gtf:
            y_start = draw_gtf_track(args, dwg, n1, n2, pos1, pos2, y_start)
        lines.append(dwg.line(
            (x_start, y_start),
            (x_start + x_scale * len(pos1), y_start),
            style='stroke:#000000;stroke-width:1',
            ))
        dwg.add(dwg.text(
            n1,
            (2*x_start + x_scale * len(pos1), y_start-10)
            ))
        dwg.add(dwg.text(
            n2,
            (2*x_start + x_scale * len(pos1), y_start+10)
            ))

        # Plot ticks
        for i in argwhere(pos1 % args.x_ticks == 1).flat:
            lines.append(dwg.line(
                (x_start + i*x_scale+.01, y_start - 0.15 * args.y_sep),
                (x_start + i*x_scale+.01, y_start + 0.15 * args.y_sep),
                style='stroke:#000000;stroke-width:1',
            ))

        if args.background:
            dwg.add(dwg.rect(
                (x_start, y_start - 5 * y_scale),
                (x_scale  * len(pos1), 10 * y_scale),
                style='stroke-width:0; fill:#BBBBBB;'
            ))
            dwg.add(dwg.line(
                (x_start, y_start - 2.5 * y_scale),
                (x_start + x_scale * len(pos1), y_start - 2.5 * y_scale),
                style='stroke:#FFFFFF;stroke-width:1;'
            ))
            dwg.add(dwg.line(
                (x_start, y_start + 2.5 * y_scale),
                (x_start + x_scale * len(pos1), y_start + 2.5 * y_scale),
                style='stroke:#FFFFFF;stroke-width:1;'
            ))
        tf_changes_i = Counter()
        for i_tf, (tf, tf_name) in enumerate(zip(args.tf, args.tf_names)):
            if args.meme_suite:
                patser1 = meme1.ix[(meme1.tf == tf) & (meme1.sequence_name == n1)]
                patser2 = meme2.ix[(meme2.tf == tf) & (meme2.sequence_name == n2)]

                patser1 = patser1.sort_values(by='score').drop_duplicates(subset='pos', keep='last')
                patser2 = patser2.sort_values(by='score').drop_duplicates(subset='pos', keep='last')

                if args.rescale_bars:
                    print(patser1.score.max())
                    top_score = max(patser1.score.max(), patser2.score.max())
                    patser1.score /= top_score / 5
                    patser2.score /= top_score / 5
                    print(patser1.score.max())

                patser1.index = patser1.pos
                patser2.index = patser2.pos

                patser1.index.name = 'pos'
                patser2.index.name = 'pos'

            else:
                in_file1 = glob(path.join(args.patser_directory, args.comp1+'*'+tf+'*'))[0]
                in_file2 = glob(path.join(args.patser_directory, args.comp2+'*'+tf+'*'))[0]

                hs1 = get_header_size(in_file1)
                if hs1:
                    patser1 = pd.read_table(in_file1, skiprows=hs1, **patser_args)
                else:
                    patser1 = pd.DataFrame(columns=patser_args['names'])
                hs2 = get_header_size(in_file2)
                if hs2:
                    patser2 = pd.read_table(in_file2, skiprows=hs2, **patser_args)
                else:
                    patser2 = pd.DataFrame(columns=patser_args['names'])

                patser1.index = [int(str(ix).strip('C')) for ix in patser1.index]
                patser2.index = [int(str(ix).strip('C')) for ix in patser2.index]
                patser1.index.name = 'pos'
                patser2.index.name = 'pos'
                patser1['pos'] = patser1.index
                patser2['pos'] = patser2.index
                patser1 = patser1.sort_values(by='pval').drop_duplicates(subset='pos')
                patser2 = patser2.sort_values(by='pval').drop_duplicates(subset='pos')
                patser1 = patser1.ix[patser1.seq == n1]
                patser2 = patser2.ix[patser2.seq == n2]
            pos = 0
            for pos in patser1.index:
                matched = has_match(pos, patser1, patser2, pos1, pos2)
                if not matched:
                    tf_changes_i['-'+tf_name] += 1
                s = float(patser1.ix[pos, 'score'])*y_scale
                r = dwg.rect(
                    (x_start + argwhere(pos1 == pos)[0,0]*x_scale, y_start-s),
                    (args.bar_width, s),
                    fill=colors[i_tf],
                    **{
                        'fill-opacity': "{:.2f}".format(1.0-match_dim*matched),
                        #'stroke-width': '0',
                        }

                    )
                r.add(svg.base.Title(tf_name + "\n" + str(patser1.ix[pos])))
                dwg_groups[tf].add(r)
            for pos in patser2.index:
                matched = has_match(pos, patser2, patser1, pos2, pos1)
                if not matched:
                    tf_changes_i['+'+tf_name] += 1
                pos = int(pos)
                s = float(patser2.ix[pos, 'score'])*y_scale
                try:
                    r = (dwg.rect(
                        (x_start + argwhere(pos2 == pos)[0,0]*x_scale, y_start),
                        (args.bar_width, s),
                        style='fill:{}; fill-opacity:{:.2f};'.format(
                            colors[i_tf],
                            1.0-match_dim*matched,
                            ),
                        ))
                    r.add(svg.base.Title(tf_name + "\n" + str(patser2.ix[pos])))
                    dwg_groups[tf].add(r)
                except IndexError:
                    pass
            pb.update(prog)
            prog += 1
        tf_changes[n1, n2] = tf_changes_i

        y_start += delta_y
    pb.finish()
    n_cols = 3
    y_start = draw_tf_cols(dwg, dwg_groups, args, y_start)
    dwg.save()
    if args.make_png:
        from subprocess import Popen
        proc = Popen(['convert', args.alignments.name+'.svg',
               args.alignments.name + '.png'])
        proc.wait()


    pprint(tf_changes)


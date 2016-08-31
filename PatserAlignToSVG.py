import svgwrite as svg
from Bio import  AlignIO, SeqIO
from argparse import ArgumentParser
from glob import glob
from os import path
from numpy import argwhere, diff, zeros
import pandas as pd
from ParseDelta import parse_records
from progressbar import ProgressBar
import itertools as it
from collections import defaultdict

seq_colors = {'T' : 'red', 'A': 'green', 'C': 'blue', 'G':'gold'}

def get_header_size(filename):
    for i, line in enumerate(open(filename)):
        if 'position=' in line:
            return i

def has_match(pos, patser1, patser2, pos1, pos2, dist=10, reltol=.2):
    score = patser1.score[pos]
    new_pos = pos2[argwhere(pos1 == pos)][0,0]
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
        bestlen, bestal = take_best(aligns_by_seq1[seq1])
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
    parser = ArgumentParser()
    parser.add_argument('--tf', '-t', default=[], nargs='*')
    parser.add_argument('--tf-names', '-T', default=[], nargs='*')
    parser.add_argument('--comp1', '-1', default='anterior')
    parser.add_argument('--comp2', '-2', default='posterior')
    parser.add_argument('--x-scale', '-x', type=float, default=.1)
    parser.add_argument('--y-scale', '-y', type=float, default=1.0)
    parser.add_argument('--y-sep', '-Y', type=float, default=50)
    parser.add_argument('--x-ticks', '-X', type=float, default=1000,
                        help='Put a marker every X-TICKS bases on the'
                        'upper-most strand')
    parser.add_argument('--bar-width', '-w', type=float, default=1)
    parser.add_argument('--show-alignment', '-A', default=False, action='store_true')
    parser.add_argument('--sequence', '-s', default=False, action='store_true',
                        help='Patser output the sequence of the match as well')
    parser.add_argument('--fasta', '-F', default=None)
    parser.add_argument('--needleall', '-N', default=False, action='store_true',
                        help='Input alignments are actually needleall output'
                        " files in EMBOSS srs format")
    parser.add_argument('alignments', type=open,
            help='File containing alignments in MUMmer .delta format')
    parser.add_argument('patser_directory', type=str)

    args =  parser.parse_args()
    if args.tf == []:
        args.tf = ['bcd', 'cad', 'gt', 'hb', 'kni', 'kr', 'hkb', 'tll', 'D', 'ftz', 'h', 'prd', 'run', 'slp1']

    if args.fasta:
        args.fasta = {rec.id: rec for rec in SeqIO.parse(args.fasta, 'fasta')}

    if args.tf_names == []:
        args.tf_names = args.tf
    elif len(args.tf_names) != len(args.tf):
        raise ValueError("Different number of TF labels and TFs: {} and {}"
                         .format(args.tf_names, args.tfs))

    return args


if __name__ == "__main__":
    args = parse_args()

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
    x_start = 10
    y_start = 1.5 * args.y_sep
    delta_y = 1.5 * args.y_sep
    match_dim = 0.5

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
                aligns_by_seq1[a1_name].append((
                    float(line.strip().split(':')[1]),
                    all_aligns[a1_name, a2_name]
                ))
        alignments = sorted(parse_best_aligns(aligns_by_seq1, max))
        n_aligns = len(alignments)
    else:
        for n_aligns, a in enumerate(parse_records(args.alignments)):
            alignments.append(a)
    pb = ProgressBar(maxval=(n_aligns+1)*len(args.tf))
    dwg = svg.Drawing(args.alignments.name+'.svg',
                     size=(x_scale * (max(len(a[0][1]) for a in alignments))
                           + 200
                           + 2 * x_start,
                           delta_y * (len(alignments)
                           * (1 + args.show_alignment) + .5)
                           + 200
                          ))
    dwg_groups = {}
    for tf in args.tf:
        dwg_groups[tf] = dwg.g(**{'class':tf})
    dwg.add(dwg.style(
        'line{stroke-width:1;stroke:#000000;}\n'
        '.hover_group{opacity:0;} \n'
        '.hover_group:hover \n'
        '{\n\topacity:1;\n'
        '\tstroke-width:1!important;'
        '\n\tstroke:#000001;\n}'
    ))

    prog = 0
    lines = []
    for ((n1, pos1), (n2, pos2)) in alignments:
        if args.show_alignment:
            fasta1 = path.join(args.patser_directory, n1 + '.fasta')
            fasta2 = path.join(args.patser_directory, n2+'.fasta')
            if path.exists(fasta1) and path.exists(fasta2):
                has_fastas = True
                seq1 = SeqIO.read(fasta1, 'fasta').seq
                seq2 = SeqIO.read(fasta2, 'fasta').seq
            elif args.fasta:
                if n1 in args.fasta and n2 in args.fasta:
                    seq1 = args.fasta[n1].seq
                    seq2 = args.fasta[n2].seq
                    has_fasta = True
                else:
                    has_fasta = False
                    seq1 = defaultdict(lambda : 'N')
                    seq2 = defaultdict(lambda : 'N')
            else:
                has_fastas = False
                seq1 = defaultdict(lambda : 'N')
                seq2 = defaultdict(lambda : 'N')

            lines.append(dwg.line(
                (x_start, y_start),
                (x_start + x_scale * len(pos1), y_start)))

            indels = diff(pos2) - diff(pos1)
            id_start = 0
            # Draw Indels
            for val, group in it.groupby(indels):
                for i, _ in enumerate(group):
                    pass
                i += 1
                dwg.add(dwg.rect(
                    (x_start + x_scale * id_start, y_start - .5 * args.y_sep * (val < 0)),
                    (x_scale * i, .5 * args.y_sep * (val != 0)),
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

        for i_tf, tf in enumerate(args.tf):
            in_file1 = glob(path.join(args.patser_directory, args.comp1+'*'+tf+'*'))[0]
            in_file2 = glob(path.join(args.patser_directory, args.comp2+'*'+tf+'*'))[0]
            patser1 = pd.read_table(in_file1, skiprows=get_header_size(in_file1), **patser_args)
            patser2 = pd.read_table(in_file2, skiprows=get_header_size(in_file2), **patser_args)
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
                s = float(patser1.ix[pos, 'score'])*y_scale
                r = dwg.rect(
                    (x_start + argwhere(pos1 == pos)[0,0]*x_scale, y_start-s),
                    (args.bar_width, s),
                    fill=colors[i_tf],
                    **{
                        'fill-opacity': "{:.2f}".format(1.0-match_dim*has_match(pos, patser1, patser2, pos1, pos2)),
                        #'stroke-width': '0',
                        }

                    )
                r.add(svg.base.Title(str(patser1.ix[pos])))
                dwg_groups[tf].add(r)
            for pos in patser2.index:
                pos = int(pos)
                s = float(patser2.ix[pos, 'score'])*y_scale
                r = (dwg.rect(
                    (x_start + argwhere(pos2 == pos)[0,0]*x_scale, y_start),
                    (args.bar_width, s),
                    style='fill:{}; fill-opacity:{:.2f};'.format(
                        colors[i_tf],
                        1.0-match_dim*has_match(pos, patser2, patser1, pos2, pos1),
                        ),
                    ))
                r.add(svg.base.Title(str(patser2.ix[pos])))
                dwg_groups[tf].add(r)
            pb.update(prog)
            prog += 1

        y_start += delta_y
    pb.finish()
    n_cols = 3
    for i, (tf, tf_name) in enumerate(zip(args.tf, args.tf_names)):
        dwg_groups[tf].add(dwg.rect(
            (x_start + 3.5 * delta_y * (i % n_cols), y_start),
            (.3*delta_y, .3*delta_y),
            **{'fill' : colors[i],
                #'stroke-width': "0",
                }
            #style='fill:{}; stroke-width:0;'.format(colors[i])

            ))
        dwg_groups[tf].add(dwg.text(tf_name,
            (x_start + 10 + 3.5 * delta_y * (i % n_cols)+.3*delta_y, y_start+0.3*delta_y)
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

    dwg.save()


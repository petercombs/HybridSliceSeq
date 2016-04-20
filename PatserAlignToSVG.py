import svgwrite as svg
from Bio import  AlignIO, SeqIO
from argparse import ArgumentParser
from glob import glob
from os import path
from numpy import zeros, argwhere, diff
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

    


def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--tf', '-t', default=[], nargs='*')
    parser.add_argument('--comp1', '-1', default='anterior')
    parser.add_argument('--comp2', '-2', default='posterior')
    parser.add_argument('--x-scale', '-x', type=float, default=.1)
    parser.add_argument('--y-scale', '-y', type=float, default=1.0)
    parser.add_argument('--bar-width', '-w', type=float, default=1)
    parser.add_argument('--show-alignment', '-A', default=False, action='store_true')
    parser.add_argument('alignments', type=open,
            help='File containing alignments in MUMmer .delta format')
    parser.add_argument('patser_directory', type=str)

    args =  parser.parse_args()
    if args.tf == []:
        args.tf = ['bcd', 'cad', 'gt', 'hb', 'kni', 'kr', 'hkb', 'tll', 'D', 'ftz', 'h', 'prd', 'run', 'slp1']

    return args


if __name__ == "__main__":
    args = parse_args()

    header_size = get_header_size(glob(path.join(args.patser_directory, '*'))[0])
    patser_args = dict(
                index_col=['pos'],
                names=['seq', '_', 'pos', '__', 'score', '___', 'pval'],
                usecols=[0,2,4,6],
                skiprows=header_size,
                na_values=[' ', '', '-', '='],
                sep=' +',
                engine='python',
                #dtype={'pos':int, 'score':float, 'pval':float},
            )

    dwg = svg.Drawing(args.alignments.name+'.svg')

    x_scale = args.x_scale
    y_scale = args.y_scale
    x_start = 10
    y_start = 20 * y_scale
    delta_y = 40 * y_scale
    match_dim = 0.9

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

    for n_aligns, a in enumerate(parse_records(args.alignments)):
        pass
    args.alignments.seek(0)
    pb = ProgressBar(maxval=(n_aligns+1)*len(args.tf))
    dwg_groups = {}
    for tf in args.tf:
        dwg_groups[tf] = dwg.g(**{'class':tf})
    dwg.add(dwg.style('line{stroke-width:1;stroke:#000000;}\n.hover_group{opacity:0;} \n.hover_group:hover \n{\n\topacity:1;\n\tstroke-width:1!important;\n\tstroke:#000000;\n}'))

    prog = 0
    lines = []
    for ((n1, pos1), (n2, pos2)) in parse_records(args.alignments):
        if args.show_alignment:
            fasta1 = path.join(args.patser_directory, n1 + '.fasta')
            fasta2 = path.join(args.patser_directory, n2+'.fasta')
            if path.exists(fasta1) and path.exists(fasta2):
                has_fastas = True
                seq1 = SeqIO.read(fasta1, 'fasta').seq
                seq2 = SeqIO.read(fasta2, 'fasta').seq
            else:
                print(path.exists(fasta1), path.exists(fasta2))
                has_fastas = False
                seq1 = defaultdict(lambda : 'N')
                seq2 = defaultdict(lambda : 'N')

            lines.append(dwg.line(
                (x_start, y_start),
                (x_start + x_scale * len(pos1), y_start)))

            indels = diff(pos2) - diff(pos1)
            id_start = 0
            for val, group in it.groupby(indels):
                for i, _ in enumerate(group):
                    pass
                i += 1
                dwg.add(dwg.rect(
                    (x_start + x_scale * id_start, y_start - 5 * y_scale * (val < 0)),
                    (x_scale * i, 5 * y_scale * (val != 0)),
                    fill="grey"
                    ))
                if val == 0:
                    for j in range(id_start, id_start + i):
                        if str(seq1[pos1[j]]) != str(seq2[pos2[j]]):
                            dwg.add(dwg.line(
                                (x_start + x_scale * j, y_start - 3 * y_scale),
                                (x_start + x_scale * j, y_start),
                                style="stroke-width:1; stroke:{};".format(seq_colors[str(seq1[pos1[j]])]),
                                ))
                            dwg.add(dwg.line(
                                (x_start + x_scale * j, y_start),
                                (x_start + x_scale * j, y_start + 3 * y_scale),
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
        for i_tf, tf in enumerate(args.tf):
            in_file1 = glob(path.join(args.patser_directory, args.comp1+'*'+tf+'*'))[0]
            in_file2 = glob(path.join(args.patser_directory, args.comp2+'*'+tf+'*'))[0]
            patser1 = pd.read_table(in_file1, **patser_args)
            patser2 = pd.read_table(in_file2, **patser_args)
            patser1 = patser1.ix[patser1.seq == n1]
            patser2 = patser2.ix[patser2.seq == n2]
            for pos in patser1.index:
                pos = int(pos)
                s = float(patser1.ix[pos, 'score'])*y_scale
                dwg_groups[tf].add(dwg.rect(
                    (x_start + argwhere(pos1 == pos)[0,0]*x_scale, y_start-s),
                    (args.bar_width, s),
                    fill=colors[i_tf],
                    **{
                        'fill-opacity': "{:.2f}".format(1.0-match_dim*has_match(pos, patser1, patser2, pos1, pos2)),
                        #'stroke-width': '0',
                        }
                    
                    ))
            for pos in patser2.index:
                pos = int(pos)
                s = float(patser2.ix[pos, 'score'])*y_scale
                dwg_groups[tf].add(dwg.rect(
                    (x_start + argwhere(pos2 == pos)[0,0]*x_scale, y_start),
                    (args.bar_width, s),
                    style='fill:{}; fill-opacity:{:.2f};'.format(
                        colors[i_tf],
                        1.0-match_dim*has_match(pos, patser2, patser1, pos2, pos1),
                        ),
                    ))
            pb.update(prog)
            prog += 1

        y_start += delta_y
    pb.finish()
    n_cols = 5
    for i, tf in enumerate(args.tf):
        dwg_groups[tf].add(dwg.rect(
            (x_start + 1.5 * delta_y * (i % n_cols), y_start),
            (.3*delta_y, .3*delta_y),
            **{'fill' : colors[i],
                #'stroke-width': "0",
                }
            #style='fill:{}; stroke-width:0;'.format(colors[i])

            ))
        dwg_groups[tf].add(dwg.text(tf,
            (x_start + 10 + 1.5 * delta_y * (i % n_cols)+.3*delta_y, y_start+0.3*delta_y)
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
            




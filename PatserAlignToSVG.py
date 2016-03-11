import svgwrite as svg
from Bio import  AlignIO
from argparse import ArgumentParser
from glob import glob
from os import path
from numpy import zeros, argwhere
import pandas as pd
from ParseDelta import parse_records
from progressbar import ProgressBar


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
            )

    dwg = svg.Drawing(args.alignments.name+'.svg')

    x_scale = args.x_scale
    y_scale = args.y_scale
    x_start = 10
    y_start = 20
    delta_y = 40 
    match_dim = 0.9

    colors=[
            'darkorange',
            'black',
            'blue',
            'gold',
            'aqua',
            'darkslateblue',
            'firebrick',
            'darkblue',
            'forestgreen',
            'lime',
            'tomato',
            'chartreuse',
            'brown',
            'crimson',
            'burlywood',
            'mediumspringgreen',
            'mediumorchid',
            'darkcyan',
            ]

    for n_aligns, a in enumerate(parse_records(args.alignments)):
        pass
    args.alignments.seek(0)
    pb = ProgressBar(maxval=(n_aligns+1)*len(args.tf))
    prog = 0
    for ((n1, pos1), (n2, pos2)) in parse_records(args.alignments):
        dwg.add(dwg.line(
            (x_start, y_start), 
            (x_start + x_scale * len(pos1), y_start),
            style='stroke:#000000;stroke-width:1',
            ))
        dwg.add(dwg.text(
            n1+'/'+n2,
            (2*x_start + x_scale * len(pos1), y_start)
            ))
        for i_tf, tf in enumerate(args.tf):
            in_file1 = glob(path.join(args.patser_directory, args.comp1+'*'+tf+'*'))[0]
            in_file2 = glob(path.join(args.patser_directory, args.comp2+'*'+tf+'*'))[0]
            patser1 = pd.read_table(in_file1, **patser_args)
            patser2 = pd.read_table(in_file2, **patser_args)
            patser1 = patser1.ix[patser1.seq == n1]
            patser2 = patser2.ix[patser2.seq == n2]
            for pos in patser1.index:
                dwg.add(dwg.rect(
                    (x_start + argwhere(pos1 == pos)[0,0]*x_scale, y_start),
                    (args.bar_width, patser1.ix[pos, 'score']*y_scale),
                    style='fill:{}; stroke-width:0; fill-opacity:{}'.format(
                        colors[i_tf],
                        1-match_dim*has_match(pos, patser1, patser2, pos1, pos2),
                        ),
                    ))
            for pos in patser2.index:
                s = patser2.ix[pos, 'score']*y_scale
                dwg.add(dwg.rect(
                    (x_start + argwhere(pos2 == pos)[0,0]*x_scale, y_start-s),
                    (args.bar_width, s),
                    style='fill:{}; stroke-width:0; fill-opacity:{}'.format(
                        colors[i_tf],
                        1-match_dim*has_match(pos, patser2, patser1, pos2, pos1),
                        ),
                    ))
            pb.update(prog)
            prog += 1

        y_start += delta_y
    pb.finish()
    n_cols = 5
    for i, tf in enumerate(args.tf):
        dwg.add(dwg.rect(
            (x_start + 100 * (i % n_cols), y_start),
            (.3*delta_y, .3*delta_y),
            style='fill:{}; stroke-width:0'.format(colors[i])
            ))
        dwg.add(dwg.text(tf,
            (x_start + 110 * (i % n_cols)+.3*delta_y, y_start+0.3*delta_y)
            ))
        y_start += delta_y * ((i%n_cols) == (n_cols-1))

    dwg.save()
            




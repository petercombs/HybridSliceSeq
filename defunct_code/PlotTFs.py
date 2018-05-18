from matplotlib.pyplot import (figure, subplot, scatter, hlines, plot, ylabel,
                               xlim, title, close, tight_layout, xticks, yticks,
                              ylim, vlines)
from numpy import ceil, log, array, isfinite
from collections import defaultdict
import pandas as pd
import fractions
from Utils import get_xs

tfs = ['cad', 'D', 'gt', 'hb', 'kni', 'prd', 'run', 'slp1']
close('all')
num_rows, num_cols = min(len(tfs), 4), int(ceil(len(tfs) / 4))
box_size = 2.5

if 'fbgn' not in locals():
    fbgns = pd.read_table('prereqs/gene_map_table_fb_2016_01.tsv', index_col=1,skiprows=5).ix[:, 0]
    gns = pd.read_table('prereqs/gene_map_table_fb_2016_01.tsv', index_col=0,skiprows=5).ix[:, 0]


ase_vals = {0.1:'11:9', 0.2: '3:2', 0.3: '13:7', 0.4: '7:3'}

def plot_ase(gene, ase, expr=None, domain_lines=None):
    xs = get_xs(ase)
    cs = ['br'[c[0]=='m'] for c in xs.index]
    print(cs)

    if gene not in ase.index:
        if gene in gns.index and gns[gene] in ase.index:
            gene = gns[gene]
        elif gene in fbgns.index and fbgns[gene] in ase.index:
            gene = fbgns[gene]
        else:
            raise KeyError("Gene {} not found!".format(gene))
    else:
        pass

    if expr is None:
        sizes = 20
    else:
        sizes = array(1.5*log(expr.ix[gene])**2)
    assert len(cs) == len(xs)
    scatter(array(xs), array(ase.ix[gene]), c=cs, s=sizes)
    yt, ytn = yticks()
    yticks(yt, [ase_val(i) for i in yt])
    ylims = ylim()
    xlims = xlim()
    hlines(0, *xlims)
    xlim(*xlims)
    if domain_lines:
        vlines(domain_lines, *ylims)
        ylim(*ylims)


def ase_val(input):
    input = abs(float(input))
    if not isfinite(input):
        return '---'
    if input == 0:
        return 'unbiased'
    A = 5 -5 * input
    B = 10 - A
    A2 = fractions.Fraction.from_float(round(100*A))
    B2 = fractions.Fraction.from_float(round(100*B))
    f = A2 / B2
    return '{}:{}'.format(f.denominator, f.numerator)

if __name__ == "__main__":
    xs = get_xs(ase)
    cs = array(['br'[c[0]=='m'] for c in xs.index])
    for gene in good_amps_peak.index.difference(good_amps_logist.index):
        figure(figsize=(1.8*box_size*num_cols, box_size*(1+num_rows)))
        for i in range(num_cols):
            subplot(num_rows + 1, num_cols, 1 + i)
            scatter(xs, ase.ix[gene], c=cs, s=1.5*log(expr.ix[gene])**2)
            hlines(0, 0, 1)
            plot(sorted(xs), peak(sorted(xs), *res_peak.ix[gene]), ':')
            ylabel(gene + '  ASE')
            xlim(0, 1)
            title(gene)
            xticks([])
            yt, ytn = yticks()
            yticks(yt, [ase_val(i) for i in yt])
            ylims = ylim()
            vlines([0.55, 0.88], *ylims)
            ylim(*ylims)

        for i, tf in enumerate(tfs):
            subplot(num_rows + 1, num_cols, 1+num_cols+i)
            scatter(xs, expr_gn.ix[tf], c = cs)
            ylabel(tf + ' FPKM')
            xlim(0, 1)
            xticks([])
            ylims = ylim()
            vlines([0.55, 0.88], *ylims)
            ylim(*ylims)
        tight_layout()


from __future__ import division
import pandas as pd
import PlotUtils as pu
from Utils import sel_startswith, pd_kwargs
from sys import argv
from os import path
from math import ceil
from matplotlib import cm

def cmap_by_prefix(prefix):

    cms = dict(
        WT = pu.ISH_CMS_5[0],
        bcd = pu.ISH_CMS_5[1],
        zld = pu.ISH_CMS_5[2],
        G20 = pu.ISH_CMS_5[3],
        hb = pu.ISH_CMS_5[4],
    )
    for p in cms:
        if prefix.startswith(p):
            return cms[p]
    return pu.ISH


if __name__ == "__main__":
    ase = pd.read_table('godot/ase_summary_by_read.tsv', **pd_kwargs)
    expr = pd.read_table('godot/summary_fb.tsv', **pd_kwargs)

    print("Read expression in")

    cdt = pd.read_table(argv[1], index_col='NAME',
                        keep_default_na=False,
                        na_values=[''])

    ase_cdt = ase.ix[cdt.index]
    exp_cdt = expr.ix[cdt.index]

    columns = (
        'melXsim',
        'simXmel',
    )

    ranges = {
        'meldominant': ('FBgn0034816', 'FBgn0250755'),
        'simdominant': ('FBgn0004087', 'FBgn0038934'),
    }

    if 'sparse' in argv[1]:
        pu.svg_heatmap(
            data=exp_cdt.select(**sel_startswith(columns)),
            filename='analysis/results/all_sparse.svg',
            norm_rows_by='max',
            progress_bar=True,
            col_sep='_sl',
            total_width=120,
            box_height=1,
            split_columns=True,
            draw_box=True,
            draw_row_labels=False,
            draw_name=True,
            cmap_by_prefix=cmap_by_prefix,
            make_hyperlinks=True,
            convert=True,
        )
        from sys import exit
        exit()
    elif '--all' in argv:
        max_height=8.5*300
        max_width=7.25*300/18
        sparsity = int(ceil(len(exp_cdt)/max_height))
        pu.svg_heatmap(
            data=exp_cdt.select(**sel_startswith(columns))[::int(sparsity)],
            filename='analysis/results/all_sparse_{:03}.svg'.format(sparsity),
            norm_rows_by='max',
            progress_bar=True,
            col_sep='_sl',
            total_width=max_width,
            box_height=1,
            split_columns=True,
            draw_box=True,
            draw_row_labels=False,
            draw_name=True,
            cmap_by_prefix=cmap_by_prefix,
            make_hyperlinks=True,
            convert=True,
        )
        pu.svg_heatmap(
            data=ase_cdt.select(**sel_startswith(columns))[::int(sparsity)],
            filename='analysis/results/all_ase_sparse_{:03}.svg'.format(sparsity),
            norm_rows_by='center0pre',
            progress_bar=True,
            col_sep='_sl',
            total_width=max_width,
            box_height=1,
            split_columns=True,
            draw_box=True,
            draw_row_labels=False,
            draw_name=True,
            cmap=cm.RdBu,
            make_hyperlinks=True,
            convert=True,
        )
        #from sys import exit
        #exit()



    for name, (gene1, gene2) in ranges.items():
        print(name)
        outname = path.join(path.dirname(argv[1]),
                            'table_{{}}_{}.svg'.format(name))
        with open(outname.replace('svg', 'txt'), 'w') as outfile:
            outfile.write('\n'.join(exp_cdt.ix[gene1:gene2].index))
            outfile.write('\n')

        pu_kwargs = dict(
            progress_bar=True,
            col_sep='_sl',
            total_width=80,
            box_height=10,
            split_columns=True,
            draw_box=True,
            draw_row_labels=True,
            draw_name=True,
            make_hyperlinks=True,
            convert=True,
        )
        pu.svg_heatmap(
            data=exp_cdt.select(**sel_startswith(columns)).ix[gene1:gene2],
            filename=outname.format('expr'),
            norm_rows_by='max',
            cmap=pu.ISH,
            **pu_kwargs
        )

        pu.svg_heatmap(
            data=ase_cdt.select(**sel_startswith(columns)).ix[gene1:gene2],
            filename=outname.format('ase'),
            norm_rows_by='center0pre',
            cmap=cm.RdBu,
            **pu_kwargs
        )


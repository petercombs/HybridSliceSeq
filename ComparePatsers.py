from Bio import AlignIO, SeqIO
from sys import argv
import pandas as pd
import numpy as np
from matplotlib.pyplot import (figure, plot, tight_layout, savefig, ylim, close, title)

def header_size(filename):
    for i, line in enumerate(open(filename)):
        if line.strip().startswith('minimum ln'):
            return i+1

if __name__ == "__main__":
    alignment = AlignIO.read(argv[1], 'emboss')
    try:
        sim = pd.read_table(argv[2], header=None,
                skiprows=header_size(argv[2]),
                sep=' ', skipinitialspace=True,
                usecols=[0,2,4,6],
                index_col=['seq', 'pos'],
                na_values=[' ', '', '-', '='],
                names=['seq', 'pos', 'score', 'pval'],
                #names=['seq', 'pos_', 'pos', 'score_', 'score', 'pval_', 'pval']
                )
        sec = pd.read_table(argv[3], header=None,
                skiprows=header_size(argv[3]),
                sep=' ', skipinitialspace=True,
                index_col=['seq', 'pos'],
                usecols=[0,2,4,6],
                na_values=[' ', '', '-', '='],
                names=['seq', 'pos', 'score', 'pval'],
                #names=['seq', 'pos_', 'pos', 'score_', 'score', 'pval_', 'pval']
                )
        patser = pd.concat([sim,sec]).dropna(how='any')
        #patser = patser.ix[patser.index.drop_duplicates()]
    except Exception as exc:
        print(argv)
        print("Failed for some reason")
        raise exc

    n = alignment.get_alignment_length()
    sim_score = np.zeros(len(alignment[0].seq.ungap('-'))+10) + np.inf
    sec_score = np.zeros(len(alignment[1].seq.ungap('-'))+10) + np.inf
    sim_pos = 0
    sec_pos = 0
    for i in range(n):
        ix_sim = (alignment[0].id, sim_pos)
        ix_sec = (alignment[1].id, sec_pos)
        if ix_sim in patser.index and ix_sec in patser.index:
            diff = patser.ix[ix_sim, 'pval'] - patser.ix[ix_sec, 'pval']
            if abs(diff) < abs(sim_score[sim_pos]):
                sim_score[sim_pos] = diff
            if abs(diff) < abs(sec_score[sec_pos]):
                sec_score[sec_pos] = diff
        sim_pos += alignment[0, i] != '-'
        sec_pos += alignment[1, i] != '-'
    if len(argv) > 4:
        sim_score.tofile(argv[4]+'_sim.txt', sep='\n')
        sec_score.tofile(argv[4]+'_sec.txt', sep='\n')
        figure(figsize=(12,2))
        plot(sim_score)
        ylim(-8,8)
        title(argv[2])
        tight_layout()
        savefig(argv[4]+'_sim.png')
        close('all')
        figure(figsize=(12,2))
        plot(sec_score)
        ylim(-8,8)
        title(argv[3])
        tight_layout()
        savefig(argv[4]+'_sec.png')
        close('all')

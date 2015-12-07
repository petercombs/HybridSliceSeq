from Bio import AlignIO, SeqIO
from sys import argv
import pandas as pd
import numpy as np
from matplotlib.pyplot import figure, plot, tight_layout, savefig, ylim

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
    sim_score = np.zeros(n)
    sec_score = np.zeros(n)
    sim_pos = 1
    sec_pos = 1
    sim_x = np.empty(n)
    sim_len = len(str(alignment[0].seq).replace('-', ''))
    for i in range(n):
        ix = (alignment[0].id, sim_pos)
        if ix in patser.index:
            sim_score[i] = patser.ix[ix, 'pval']
        else:
            sim_score[i] = np.nan
        ix = (alignment[1].id, sec_pos)
        if ix in patser.index:
            sec_score[i] = patser.ix[ix, 'pval']
        else:
            sec_score[i] = np.nan
        sim_pos += alignment[0, i] != '-'
        sec_pos += alignment[1, i] != '-'
        sim_x[i] = sim_pos
    if len(argv) > 4:
        (sim_score-sec_score)[1000:].tofile(argv[4]+'.txt', sep='\n')
        figure(figsize=(12,2))
        plot(sim_x[1000:], (sim_score-sec_score)[1000:])
        ylim(-8,8)
        tight_layout()
        savefig(argv[4]+'.png')

        



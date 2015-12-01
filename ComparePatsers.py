from Bio import AlignIO, SeqIO
from sys import argv
import pandas as pd
import numpy as np

if __name__ == "__main__":
    alignment = AlignIO.read(argv[1], 'emboss')
    sim = pd.read_table(argv[2], header=None, skiprows=39, 
            sep=' ', skipinitialspace=True, 
            names=['seq', 'pos_', 'pos', 'score_', 'score', 'pval_', 'pval'])
    sec = pd.read_table(argv[3], header=None, skiprows=39, 
            sep=' ', skipinitialspace=True, 
            names=['seq', 'pos_', 'pos', 'score_', 'score', 'pval_', 'pval'])
    patser = pd.concat([sim,sec], ignore_index=True)

    n = alignment.get_alignment_length()
    sim_score = np.zeros(n)
    sec_score = np.zeros(n)
    sim_pos = 1
    sec_pos = 1
    for i in range(n):
        q = patser.query('seq == "{}" and pos == {}' .format(alignment[0].id, sim_pos))
        sim_score[i] = q.score if len(q) else np.nan
        q = patser.query('seq == "{}" and pos == {}' .format(alignment[1].id, sec_pos))
        sec_score[i] = q.score if len(q) else np.nan
        sim_pos += alignment[0, i] != '-'
        sec_pos += alignment[1, i] != '-'



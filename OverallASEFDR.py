import numpy as np
import pandas as pd
from progressbar import ProgressBar as pbar
from Utils import (sel_startswith, pd_kwargs, get_chroms)
from multiprocessing import Pool
from GetASEStats import FRAC_FOR_MATERNAL, EXPR_MIN

def get_randomized_scores(ase):
    rand = np.random.randint(2, size=ase_expr.shape)*2-1
    flipped = ase * rand
    melXsim = flipped.select(**sel_startswith('melXsim'))
    simXmel = flipped.select(**sel_startswith('simXmel'))
    weaker_sim_bias = np.min([melXsim.T.quantile(1-FRAC_FOR_MATERNAL),
                              simXmel.T.quantile(1-FRAC_FOR_MATERNAL)],
                             axis=0)
    weaker_mel_bias = np.max([melXsim.T.quantile(FRAC_FOR_MATERNAL),
                              simXmel.T.quantile(FRAC_FOR_MATERNAL)],
                             axis=0)
    return (weaker_mel_bias, weaker_sim_bias)


if __name__ == "__main__":
    ase = locals().get('ase', None)
    expr = locals().get('expr', None)
    if ase is None or expr is None or not np.all(ase.index == expr.index):
        print("reloading files")
        expr = pd.read_table('analysis_godot/summary.tsv', **pd_kwargs).dropna(how='all', axis=1)
        ase = (pd
               .read_table('analysis_godot/ase_summary_by_read.tsv',
                           **pd_kwargs
                           )
               .dropna(how='all', axis=1)
               .dropna(how='all', axis=0)
               .select(**sel_startswith(('melXsim', 'simXmel')))
              )
        chrom_of = get_chroms()

        males = ('melXsim_cyc14C_rep3', 'simXmel_cyc14C_rep2')
        on_x = [chrom_of[gene] == 'X' if gene in chrom_of else False for gene in ase.index]
        is_male = [col.startswith(males) for col in ase.columns]
        ase.ix[on_x, is_male] = np.nan

    melXsim_expr = expr.select(**sel_startswith('melXsim'))
    simXmel_expr = expr.select(**sel_startswith('simXmel'))
    melXsim_ase = ase.select(**sel_startswith('melXsim'))
    simXmel_ase = ase.select(**sel_startswith('simXmel'))
    melXsim_is_expr = (melXsim_expr > EXPR_MIN)
    simXmel_is_expr = (simXmel_expr > EXPR_MIN)
    all_is_expr = expr > EXPR_MIN

    min_per_crossdir = 10
    expr_both = (
        (melXsim_ase.T.count() > min_per_crossdir)
        & (simXmel_ase.T.count() > min_per_crossdir)
        & (melXsim_is_expr.T.sum() > min_per_crossdir)
        & (simXmel_is_expr.T.sum() > min_per_crossdir)
    )
    ase_expr = ase.ix[expr_both]
    print("Found {} good genes".format(len(ase_expr)))
    n_reps = 1000
    mel_biases = pd.DataFrame(index=ase_expr.index, columns=range(n_reps))
    sim_biases = pd.DataFrame(index=ase_expr.index, columns=range(n_reps))
    with Pool() as p:
        results = [None for i in range(n_reps)]
        for i in range(n_reps):
            results[i] = p.apply_async(get_randomized_scores, (ase_expr, ))

        for i, res in pbar(max_value=n_reps)(enumerate(results)):
            mel_bias, sim_bias = res.get()
            mel_biases.ix[:, i] = mel_bias
            sim_biases.ix[:, i] = sim_bias





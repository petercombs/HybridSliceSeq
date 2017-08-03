from __future__ import print_function
import pandas as pd
import numpy as np
import Utils as ut
import sys
from scipy.stats.distributions import chi2
from bisect import bisect
from progressbar import ProgressBar as pb
from multiprocessing import Pool


def ase_score_to_prob(ase_score):
    return np.log2(2/(1-ase_score) - 1)



def fishers_method(pvals, clipped=False):
    chi2_2k = -2 * np.sum(np.log(np.array(pvals)))
    return chi2.sf(chi2_2k, 2*len(pvals))

def shuffle_ase(good_ase, embs):
    ase_shuffle = good_ase * (np.random.randint(2, size=good_ase.shape)*2-1)
    return {
        emb: (ase_shuffle.T
              .select(ut.startswith(emb))
              .sum()
             )
        for emb in embs
    }



MIN_SLICES_FOR_ASE = 50
N_SHUFFLES = int(5e4)

if __name__ == "__main__":
    ase = (pd.read_table('analysis_godot/ase_summary_by_read.tsv', **ut.pd_kwargs)
           .select(**ut.sel_startswith(('melXsim', 'simXmel')))
           .drop('---', axis=1, errors='ignore')
          )


    good_ase = (ase
                .copy()
                #.applymap(ase_score_to_prob)
                .ix[ase.T.count() > MIN_SLICES_FOR_ASE]
               )
    if 'syns' not in locals():
        syns = ut.get_synonyms()
        chrom_of = ut.get_chroms(syns)
        males = ('melXsim_cyc14C_rep3', 'simXmel_cyc14C_rep2')
        on_x = pd.Series({gene: chrom_of[gene] == 'X'
                          for gene in good_ase.index})
        is_male = [col.startswith(males) for col in good_ase.columns]
    good_ase.ix[on_x, is_male] = pd.np.nan

    embs = {col.split('_sl')[0] for col in ase.columns}
    by_emb_sums = pd.DataFrame(index=good_ase.index, columns=sorted(embs))
    for emb in embs:
        by_emb_sums.ix[:, emb] = good_ase.T.select(ut.startswith(emb)).sum()

    rando_sums = {
        emb: pd.DataFrame(index=good_ase.index, columns=np.arange(N_SHUFFLES),
                          data=np.nan)
        for emb in embs
    }
    print("Randomizing", file=sys.stderr)
    with Pool() as p:
        results = []
        for i in range(N_SHUFFLES):
            results.append(p.apply_async(shuffle_ase, (good_ase, embs)))

        for i in pb()(list(range(N_SHUFFLES))):
            result = results[i].get()
            for emb in result:
                rando_sums[emb].ix[:, i] = result[emb]

    pvals_by_emb = pd.DataFrame(index=good_ase.index, columns=sorted(embs),
                                dtype=float)
    for gene in pb()(pvals_by_emb.index):
        for emb in pvals_by_emb.columns:
            pvals_by_emb.ix[gene, emb] = (
                (
                    bisect(
                        sorted(rando_sums[emb].ix[gene]),
                        by_emb_sums[emb].ix[gene]
                    )
                )
                / N_SHUFFLES
            )
    pvals_by_emb = pvals_by_emb.applymap(
        lambda x: 1-x if 1-x < x else x
    )
    pvals_by_emb = pvals_by_emb.applymap(
        lambda x: x if x else 1 / (2*N_SHUFFLES)
    )


    combined_pvals = pd.Series(index=good_ase.index)
    for gene in combined_pvals.index:
        if on_x[gene]:
            combined_pvals.ix[gene] = fishers_method(pvals_by_emb.ix[gene].drop(list(males)))
        else:
            combined_pvals.ix[gene] = fishers_method(pvals_by_emb.ix[gene])

    # Note the double negation to fix instances of NaN's, a common source of
    # which is X genes in males
    inconsistent_dir = ~(
        (~(np.sign(by_emb_sums).T != 1.0)).all()
        | (~(np.sign(by_emb_sums).T != -1.0)).all()
    )

    combined_pvals.ix[inconsistent_dir] = 1

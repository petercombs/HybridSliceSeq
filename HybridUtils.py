from matplotlib.pyplot import figure
import numpy as np
import pandas as pd
from Utils import startswith
try:
    from CompareDESeqs import ase_bars2
    from Annote import AnnoteFinder
    ab2 = True
except ImportError:
    ab2 = False

from scipy.stats import ttest_1samp

def get_classes(ase, **kwargs):
    pbar = kwargs.pop('pbar', lambda : lambda x: x)
    classes = pd.DataFrame(index=ase.index,
                           columns={col.split('_')[0] for col in ase.columns},
                           data=np.nan)
    for gene in pbar()(ase.index):
        classes.ix[gene] = is_directionally_biased(ase, gene, **kwargs)
    return classes

def is_maternal(ase, gene, **kwargs):
    mat_direction = np.array([1 if col.startswith('sim') else -1
                     for col in ase.columns])
    biases =  is_directionally_biased(ase, gene, mat_direction,
                                      too_few_slices_val=0, **kwargs)
    return np.all(list(k==1 for k in biases.values()))

def is_paternal(ase, gene, **kwargs):
    pat_direction = np.array([-1 if col.startswith('sim') else 1
                     for col in ase.columns])
    biases = is_directionally_biased(ase, gene, pat_direction, **kwargs)
    return np.all(list(biases.values()))

def is_sim(ase, gene, **kwargs):
    sim_direction = np.array([1 for col in ase.columns])
    biases = is_directionally_biased(ase, gene, sim_direction, **kwargs)
    return np.all(list(biases.values()))

def is_zyg(ase, gene, **kwargs):
    biases = is_directionally_biased(ase, gene **kwargs)
    return not np.any(list(biases.values()))



def is_directionally_biased(ase, gene, bias_direction=None, style='ttest', ase_level=0.33,
                            min_slices=10, too_few_slices_val=99,
                            frac_for_biased=0.65, two_tailed=False, alpha=.05):
    if bias_direction is None:
        bias_direction = [1 for col in ase.columns]
    genotypes = {col.split('_')[0] for col in ase.columns}
    biases = {}
    for genotype in genotypes:
        genease = (ase.ix[gene] * bias_direction).select(startswith(genotype))
        if style == 'ttest':
            tstat, pval = ttest_1samp(genease, 0, nan_policy='omit')
            if isinstance(pval, np.ma.core.MaskedConstant):
                biases[genotype] = too_few_slices_val
                continue
            if two_tailed:
                biases[genotype] = np.sign(tstat) * (pval * len(ase) < alpha)

            else:
                pval = pval/2 if tstat > 0 else 1-pval/2
                biases[genotype] = pval * len(ase)  < alpha

        elif style == 'cutoff':
            slices_with_aseval = genease.count()
            if slices_with_aseval < min_slices:
                biases[genotype] = too_few_slices_val
                continue
            biases[genotype] = 0
            for dir in [-1, 1]:
                if ((dir * genease > ase_level).sum()
                    > max(frac_for_biased * slices_with_aseval, min_slices)):
                    biases[genotype] = dir
                    break
        else:
            raise NotImplementedError("Don't know how to use test style '{}'".format(style))
    return biases




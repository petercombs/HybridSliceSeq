import DistributionDifference as dd
import pandas as pd
import Utils as ut
import itertools as it
from progressbar import ProgressBar as pb
from builtins import sum


EXPR_MIN = 5

if __name__ == "__main__":
    expr = (pd
            .read_table('analysis_godot/summary.tsv', **ut.pd_kwargs)
            .drop('---', axis=1)
           )

    mel = expr.select(**ut.sel_startswith('melXmel_'))
    sim = expr.select(**ut.sel_startswith('simXsim_'))
    hyb = expr.select(**ut.sel_startswith(('melXsim', 'simXmel')))
    expr_in_mel = (mel.max(axis=1) > EXPR_MIN)
    expr_in_sim = sim.max(axis=1) > EXPR_MIN
    expr_in_hybrids = (hyb.max(axis=1) > EXPR_MIN)
    expr_in_all = (expr_in_mel & expr_in_sim & expr_in_hybrids)

    expr = expr.ix[expr_in_all]

    embryo_types = {c.split('_sl')[0].split('_rep')[0] for c in expr.columns}
    embryos = {}
    for etype in embryo_types:
        embryos[etype] = {c.split('_sl')[0]
                          for c in expr.columns if
                          c.startswith(etype)}

    combs = sum(
        [sorted(it.combinations(e,2)) for e in embryos.values()],
        []
    )
    combs += list(it.product(embryos['melXsim_cyc14C'],
                             embryos['simXmel_cyc14C']))
    emds = pd.DataFrame(
        index=expr.index,
        columns=["{}-{}".format(*c) for c in combs],
        data=-1
    )
    for gene in pb()(expr.index):
        for e1, e2 in combs:
            emds.ix[gene, "{}-{}".format(e1, e2)] = (
                dd.earth_mover_multi_rep(
                    expr.ix[gene].select(ut.startswith(e1))+EXPR_MIN,
                    expr.ix[gene].select(ut.startswith(e2))+EXPR_MIN,
                ))



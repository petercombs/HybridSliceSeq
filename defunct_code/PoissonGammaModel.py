import pymc3 as pm
import pandas as pd
import Utils as ut

def is_mat(label):
    mother = label.split('X')[0]
    refalt = label.split('_')[-1]
    return ((mother == 'mel' and refalt == 'ref')
            or (mother == 'sim' and refalt == 'alt'))

endswith = lambda y: lambda x: x.endswith(y)

if __name__ == "__main__":
    refalt_data = (pd.read_table('analysis_godot/ase_summary_refalt.tsv',
                                **ut.pd_kwargs)
                   .select(**ut.sel_startswith(('melXsim', 'simXmel')))
                  )

    samples = [i.rsplit('_', 1)[0] for i in refalt_data.columns]
    for gene in ['skpA']:
        with pm.Model() as model:
            mu = pm.Gamma('mu', 1/2, 1/2)
            alpha = pm.Gamma('alpha', 1/2, 1/2)
            q = 0.5
            n = 0
            for n, i in enumerate(samples):
                beta_i = pm.Gamma('beta_{}'.format(n), 1/2, 1/2)
                x_i = pm.Poisson('alt_{}'.format(n), mu * alpha * beta_i * q,
                                 observed=refalt_data.ix[gene]
                                 .select(endswith('alt'))
                                 .select(ut.startswith(i))
                                )
                y_i = pm.Poisson('ref_{}'.format(n), mu * alpha * beta_i * q,
                                 observed=refalt_data.ix[gene]
                                 .select(endswith('ref'))
                                 .select(ut.startswith(i))
                                )

            print("Finding MAP")
            start = pm.find_MAP(progressbar=True)
            step = pm.Metropolis()
            trace = pm.sample(20000, step, start=start, progressbar=True)






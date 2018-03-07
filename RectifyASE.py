import Utils as ut
import pandas as pd


if __name__ == "__main__":
    ase = (pd.read_table('analysis_godot/ase_summary_by_read.tsv', **ut.pd_kwargs)
           .dropna(how='all', axis=1)
           .select(**ut.sel_startswith(('melXsim', 'simXmel')))
          )

    syns = ut.get_synonyms()
    chrom_of = ut.get_chroms(syns)

    males = ('melXsim_cyc14C_rep3', 'simXmel_cyc14C_rep2')
    on_x = [chrom_of[gene] == 'X' for gene in ase.index]
    is_male = [col.startswith(males) for col in ase.columns]
    ase_nomaleX = ase.copy()
    ase_nomaleX.ix[on_x, is_male] = pd.np.nan
    ase = ase_nomaleX

    rectified_ase = ase.multiply([1 if ix.startswith('melXsim') else -1
                                  for ix in ase.columns])


    rectified_ase.to_csv('analysis_godot/rectified_ase_summary_by_read.tsv',
                        sep='\t', na_rep='---', float_format='%.04f')


import pandas as pd


if __name__ == "__main__":
    if 'ase' not in locals() or ('reload_ase' in locals() and locals()['reload_ase']):
        print("Reloading data")

        ase = pd.read_table('analysis_godot/ase_summary.tsv', **kwargs).dropna(how='all', axis=1)
        expr = pd.read_table('analysis_godot/summary_fb.tsv', **kwargs).dropna(how='all', axis=1)

        lott = pd.read_table('prereqs/journal.pbio.1000590.s002', index_col=0, keep_default_na=False, na_values=[''])
        lott['fbgn'] = to_fbgn[lott.index]
        lott.drop_duplicates(subset='fbgn', keep=False, inplace=True)
        lott.index = lott.fbgn

    lott_mat_only = None

from os import path
from numpy import shape, linspace, sum, isfinite
import pandas as pd


def get_bam_length(samfile):
    start = samfile.tell()
    maxval = path.getsize(samfile.filename) * 2**16
    # I don't know why that 2**16 factor is there!
    return maxval + 2**16, start


def strip_to_number(dataval, chars='\'" \t #'):
    return to_number(dataval.strip(chars))


def to_number(dataval):
    """ A forgiving number converter.

    Will convert to int if possible, float otherwise, and if neither, will
    return the input.
    """
    try:
        datavalf = float(dataval)
        # If we could convert it to a float, it might have been an
        # int
        try:
            return int(dataval)
        except ValueError:
            # not an int, but since we got to the inner try, it is a
            # float
            return datavalf

    except ValueError:
        return dataval


def contains(string_or_iterable):
    if isinstance(string_or_iterable, str):
        return lambda x: string_or_iterable in x
    else:
        return lambda x: any(i in x for i in string_or_iterable)

def startswith(string_or_iterable):
    if not isinstance(string_or_iterable, str):
        string_or_iterable = tuple(string_or_iterable)
    return lambda x: x.startswith(string_or_iterable)

def sel_contains(string_or_iterable):
    return dict(crit=contains(string_or_iterable), axis=1)

def sel_startswith(string_or_iterable):
    return dict(crit=startswith(string_or_iterable), axis=1)

def center_of_mass(data):
    if 'columns' in dir(data) and 'rep' in data.columns[0]:
        reps = {c.split('_sl')[0] for c in data.columns}
        retval = 0
        for rep in reps:
            retval += center_of_mass_onerep(data.select(**sel_startswith(rep)))
        return retval / len(reps)
    elif 'index' in dir(data) and 'rep' in data.index[0]:
        reps = {c.split('_sl')[0] for c in data.index}
        retval = 0
        for rep in reps:
            retval += center_of_mass_onerep(data.select(startswith(rep)))
        return retval / len(reps)

    else:
        return center_of_mass_onerep(data)

def center_of_mass_onerep(data):
    dims = shape(data)
    cols = dims[-1]
    xs = linspace(0, 1, cols, endpoint=True)
    data_clean = data.copy()
    data_clean[~isfinite(data_clean)] = 0
    data_clean += 0.01
    return sum(data_clean * xs, axis=len(dims)-1)/sum(data_clean, axis=len(dims)-1)


def load_to_locals(locals, expr_min=15, fname='analysis/summary.tsv'):
    read_table_args = dict(keep_default_na=False, na_values=['---', ''], index_col=0)

    # These are manually, and empirically determined.
    bad_cols = (
        'bcd_cyc14D_rep2_sl06_FPKM',
        'bcd_cyc14D_rep2_sl16_FPKM',
        'bcd_cyc14D_rep1_sl14_FPKM',
        #'WT_cyc14D_sl15_FPKM',
        'G20_cyc14D_rep1_sl08_FPKM',
    )

    all_expr = (pd.read_table(fname, **read_table_args)
                .sort_index())
    for col in bad_cols:
        if col in all_expr.columns:
            all_expr.ix[:, col] = pd.np.nan
        else:
            print("Column {} is missing")
    top_expr = all_expr.max(axis=1)
    all_expr = all_expr.ix[top_expr > expr_min]
    wt  = all_expr.select(**sel_startswith('WT'))
    bcd = all_expr.select(**sel_startswith('bcd'))
    zld = all_expr.select(**sel_startswith('zld'))
    g20 = all_expr.select(**sel_startswith('G20'))
    hb  = all_expr.select(**sel_startswith('hb'))
    locals['all_expr'] = all_expr
    locals['wt'] = wt
    locals['bcd'] = bcd
    locals['g20'] = g20
    locals['zld'] = zld
    locals['hb'] = hb

    by_cycle = {}
    for sub_df_name in 'wt bcd zld g20 hb'.split():
        sub_df = locals[sub_df_name]
        cycs = {col.split('_sl')[0].split('_',1)[1] for col in sub_df.columns}
        cycs.update({col.split('_')[1] for col in sub_df.columns})
        cyc_embs = {}
        by_cycle[sub_df_name] = cyc_embs
        for cyc in cycs:
            cyc_embs[cyc] = sub_df.select(**sel_contains(cyc))
        locals[sub_df_name+'s'] = cyc_embs
    return (all_expr,
            [wt, bcd, zld, g20, hb, ],
            [locals[i] for i in 'wts bcds zlds g20s hbs'.split()],
            by_cycle)

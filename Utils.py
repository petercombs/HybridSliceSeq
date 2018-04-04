import sys
from os import path
from numpy import shape, linspace, sum, isfinite
import pandas as pd
from collections import defaultdict
from glob import glob


def get_bam_length(samfile):
    start = samfile.tell()
    maxval = path.getsize(samfile.filename) * 2**16
    # I don't know why that 2**16 factor is there!
    return maxval + 2**16, start


def strip_to_number(dataval, chars='\'" \t #'):
    return to_number(dataval.strip(chars))

def true_index(series):
    'Returns elements in the index of the series where the element is true'
    return series.index[series.astype(bool)]


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

def not_contains(string_or_iterable):
    return lambda x: not contains(string_or_iterable)(x)

def startswith(string_or_iterable):
    if not isinstance(string_or_iterable, str):
        string_or_iterable = tuple(string_or_iterable)
    return lambda x: x.startswith(string_or_iterable)

def not_startswith(string_or_iterable):
    return lambda x: not startswith(string_or_iterable)(x)

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

for dir_name in sys.path:
    fname = path.join(dir_name, 'prereqs/gene_map_table_fb_2014_03.tsv')
    if path.exists(fname):
        fbgns = pd.read_table(fname,
                              keep_default_na=False, na_values=['---'],
                              skipfooter=2, engine='python',
                              index_col=1,skiprows=5).ix[:, 0]
        gn_to_fbgn = {val: key for key, val in fbgns.items()}
        break
else:
    pass
    #raise ImportError('Could not find Gene mapping table')

def get_synonyms():
    gn_to_fbgn = defaultdict(lambda : 'NOTPRESENT')
    file = [path.join(dirname, fname)
            for fname in ['fbgn_synonyms.tsv', 'synonyms.tsv']
            for dirname in ['.', 'Reference', 'prereqs']
            if path.exists(path.join(dirname, fname))
           ][0]
    secondaries = {}
    for line in open(file):
        line = line.strip().split('\t')
        if len(line) < 2: continue
        fbgn = line[0]
        gn_to_fbgn[line[1]] = fbgn
        if len(line) == 3 or len(line) == 4:
            continue
        for synonym in line[-1].split(','):
            secondaries[synonym] = fbgn
    for synonym in secondaries:
        if synonym in gn_to_fbgn:
            # Don't clobber primary names with the synonym
            continue
        gn_to_fbgn[synonym] = secondaries[synonym]
    gn_to_fbgn['Dec1'] = gn_to_fbgn['dec-1']
    return pd.Series(gn_to_fbgn)

def load_to_locals(locals, expr_min=15):
    read_table_args = dict(keep_default_na=False, na_values=['---', ''], index_col=0)

    # These are manually, and empirically determined.
    bad_cols = (
        'bcd_cyc14D_rep2_sl06_FPKM',
        'bcd_cyc14D_rep2_sl16_FPKM',
        'bcd_cyc14D_rep1_sl14_FPKM',
        'WT_cyc14D_sl15_FPKM',
        'G20_cyc14D_rep1_sl08_FPKM',
    )

    all_expr = (pd.read_table('analysis/summary.tsv', **read_table_args)
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


pd_kwargs = dict(
    index_col=0, keep_default_na=False, na_values=['---', '', '-'],
)

def get_xs(dataframe):
    max_slice = defaultdict(int)
    if isinstance(dataframe, pd.Series):
        dataframe = pd.DataFrame(dataframe).T
    for sl in dataframe.columns:
        sl = sl.split('_sl')
        emb = sl[0]
        if emb == '---':
            max_slice[emb] = 0
        else:
            max_slice[emb] = max(max_slice[emb], int(sl[1][0:2]))

    return pd.Series(index=dataframe.select(**sel_contains('sl')).columns,
                   data=[(int(a.split('_sl')[1][:2])-1)/(max_slice[a.split('_sl')[0]]-1)
                         for a in dataframe.columns if 'sl' in a])

def get_nearest_slice(query, other_frame):
    other_xs = get_xs(other_frame)
    return abs(query - other_xs).idxmin()

def get_coords(syns=None):
    coord_of = {}
    for fname in reversed(sorted(glob('prereqs/gene_map_table*.tsv'))):
        for row in open(fname):
            if row.startswith('#') or not row.strip():
                continue
            data = row.split()
            try:
                coord = int(data[-1].split(':')[1].split('..')[0])
                if data[1] not in coord_of:
                    coord_of[data[1]] = coord
                if data[0] not in coord_of:
                    coord_of[data[0]] = coord
            except:
                continue
    if syns is not None:
        for key in syns.index:
            if key not in coord_of and syns[key] in coord_of:
                coord_of[key] = coord_of[syns[key]]
    return pd.Series(coord_of)


def get_chroms(syns=None):
    chrom_of = {}
    for fname in reversed(sorted(glob('prereqs/gene_map_table*.tsv'))):
        for row in open(fname):
            if row.startswith('#') or not row.strip():
                continue
            data = row.split()
            chrom = data[-1].split(':')[0]
            if data[1] not in chrom_of:
                chrom_of[data[1]] = chrom
            if data[0] not in chrom_of:
                chrom_of[data[0]] = chrom
    if syns is not None:
        for key in syns.index:
            if key not in chrom_of and syns[key] in chrom_of:
                chrom_of[key] = chrom_of[syns[key]]
    return pd.Series(chrom_of)

def parse_annotation(gtf_annote_field):
    recs = gtf_annote_field.strip().strip(';').split(';')
    retval = {}
    for rec in recs:
        key, val = rec.strip().split(' ', 1)
        retval[key] = val.strip('"\'; ')
    return retval


def get_exons_before_CDS(exon_parts, ref_annotation,
                         progress=False):
    epn = 'exonic_part_number'
    if progress:
        from tqdm import tqdm
    CDS_starts = {}
    for line in tqdm(open(ref_annotation)):
        data = line.strip().split('\t')
        if data[2] != 'CDS': continue
        annots = parse_annotation(data[-1])
        if data[6] == '+':
            CDS_starts[annots['gene_id']] = min(int(data[3]),
                                                CDS_starts.get('gene_id', 1e10))
        elif data[6] == '-':
            CDS_starts[annots['gene_id']] = max(int(data[4]),
                                                CDS_starts.get('gene_id', 0))
            pass
        else:
            raise ValueError('Strand field has unknown type "{}"'
                             .format(data[6]))
    pre_exons = defaultdict(list)
    strands = {}
    for line in tqdm(open(exon_parts)):
        data = line.strip().split('\t')
        if data[2] != 'exonic_part': continue
        annots = parse_annotation(data[-1])
        for gene in annots['gene_id'].split('+'):
            strand = data[6]
            start, stop = int(data[3]), int(data[4])
            strands[gene] = strand
            if strand == '+':
                if start < CDS_starts[gene]:
                    pre_exons[gene].insert(0, annots[epn])
            elif strand == '-':
                if stop > CDS_starts[gene]:
                    pre_exons[gene].append(annots[epn])
    return pre_exons




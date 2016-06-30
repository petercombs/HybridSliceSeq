from matplotlib.pyplot import figure
import pandas as pd
try:
    from CompareDESeqs import ase_bars2
    from Annote import AnnoteFinder
    ab2 = True
except ImportError:
    ab2 = False


startswith = lambda y: lambda x: x.startswith(y)

if 'expr' not in locals():
    expr = pd.read_table('summary.tsv', index_col=0)

if 'go' not in locals():
    go = pd.read_table('gene_association.fb', header=None, skiprows=5,
                       names=['FB', 'FBgn', 'symbol', 'blank', 'go', 'ref', 'type', 'blank2', 'valid',
                              'longname', 'huh', 'product', 'taxon', 'number', 'source', 'blank3',
                              'blank4'],
                       sep='\t', comment='!',
                       index_col=None)
if 'orthologs' not in locals():
    orthologs = pd.read_table('gene_orthologs_fb_2015_03.tsv',
                              index_col='Ortholog_FBgn_ID', header=3,
                              names=['FBgn_ID', 'GeneSymbol',
                                     'Arm/Scaffold', 'Location', 'Strand',
                                     'Ortholog_FBgn_ID',
                                     'Ortholog_GeneSymbol',
                                     'Ortholog_Arm/Scaffold',
                                     'Ortholog_Location', 'Ortholog_Strand',
                                     'OrthoDB_Group_ID'])
    gnix_orthologs = pd.read_table('gene_orthologs_fb_2015_03.tsv',
                                   index_col='Ortholog_GeneSymbol', header=3,
                                   names=['FBgn_ID', 'GeneSymbol',
                                          'Arm/Scaffold', 'Location', 'Strand',
                                          'Ortholog_FBgn_ID',
                                          'Ortholog_GeneSymbol',
                                          'Ortholog_Arm/Scaffold',
                                          'Ortholog_Location', 'Ortholog_Strand',
                                          'OrthoDB_Group_ID'])
    melix_orthologs = pd.read_table('gene_orthologs_fb_2015_03.tsv',
                              index_col='FBgn_ID', header=3,
                              names=['FBgn_ID', 'GeneSymbol',
                                     'Arm/Scaffold', 'Location', 'Strand',
                                     'Ortholog_FBgn_ID',
                                     'Ortholog_GeneSymbol',
                                     'Ortholog_Arm/Scaffold',
                                     'Ortholog_Location', 'Ortholog_Strand',
                                     'OrthoDB_Group_ID'])
    melix_orthologs =(
        melix_orthologs.ix[melix_orthologs
                           .Ortholog_GeneSymbol
                           .apply(startswith('Dsim'))]
    )
    melgnix_orthologs = pd.read_table('gene_orthologs_fb_2015_03.tsv',
                              index_col='GeneSymbol', header=3,
                              names=['FBgn_ID', 'GeneSymbol',
                                     'Arm/Scaffold', 'Location', 'Strand',
                                     'Ortholog_FBgn_ID',
                                     'Ortholog_GeneSymbol',
                                     'Ortholog_Arm/Scaffold',
                                     'Ortholog_Location', 'Ortholog_Strand',
                                     'OrthoDB_Group_ID'])
    melgnix_orthologs =(
        melgnix_orthologs.ix[melgnix_orthologs
                           .Ortholog_GeneSymbol
                           .apply(startswith('Dsim'))]
    )

if 'adj_ases' not in locals():
    chroms = pd.read_pickle('results/chroms.pkl')
    tmp = pd.read_table('oefemale_gene_ase.tsv', index_col=0)
    ases = pd.DataFrame(index=tmp.index)
    for tissue in 'oe fb'.split():
        for sex in ['male', 'female']:
            for rep in ['', '_r2']:
                fname = tissue+sex+rep+'_counts_10_nomin_gene_ase.tsv'
                tmp = pd.read_table(fname, index_col=0)
                ases[tissue+sex+rep+"_REF"] = tmp.REFERENCE_COUNTS
                ases[tissue+sex+rep+"_ALT"] = tmp.ALT_COUNTS
    simchr = chroms.ix[ases.index]
    notx = ((simchr != 'Scf_X')
            & (simchr != 'Scf_NODE_24415')
            & (simchr != 'Scf_NODE_39788')
            & (simchr != 'Scf_NODE_111243')
            &(simchr != 'Scf_NODE_12414')
            &(simchr != 'Scf_NODE_13260')
            &(simchr != 'Scf_NODE_13320')
            &(simchr != 'Scf_NODE_14159')
            &(simchr != 'Scf_NODE_15140')
            &(simchr != 'Scf_NODE_19857')
            &(simchr != 'Scf_NODE_22823')
            &(simchr != 'Scf_NODE_24415')
            &(simchr != 'Scf_NODE_25881')
            &(simchr != 'Scf_NODE_29806')
            &(simchr != 'Scf_NODE_3296')
            &(simchr != 'Scf_NODE_3417')
            &(simchr != 'Scf_NODE_37002')
            &(simchr != 'Scf_NODE_39788')
            &(simchr != 'Scf_NODE_53243')
            &(simchr != 'Scf_NODE_62956')
            &(simchr != 'Scf_NODE_68540')
            &(simchr != 'Scf_NODE_75619'))
    libsizes_dedup = pd.Series(index=ases.columns)
    libsizes = pd.Series.from_csv('libsizes.tsv', sep='\t')
    for row in libsizes_dedup.index:
        libsizes_dedup.ix[row] = libsizes.ix['analysis/on_sim/{}_STAR_RNASEQ_wasp_dedup.bam'.format(row.rsplit('_', 1)[0])]
    adj_ases = ( ases
                .ix[notx]
                .divide(libsizes_dedup)
                .replace(to_replace=pd.np.nan, value=0)
                * 1e6 + 1
               )

sim_orthologs = (orthologs.ix[adj_ases.index]
                 .dropna(how='all')
                 .drop_duplicates(subset='FBgn_ID'))

sim_orthologs['Ortholog_loc_lo'] = [int(loc.split('..')[0]) for loc in sim_orthologs.Ortholog_Location]



if ab2:
    ab2 = lambda x: [figure(figsize=(3,2)),
                     ase_bars2(x, adj_ases+1, 'oefemale', False,
                               sim_orthologs.drop_duplicates(subset='Ortholog_GeneSymbol').GeneSymbol.get(x, 'No Ortholog'),
                               expr)]
else:
    ab2 = lambda x: False


if 'map_info' not in locals():
    import re
    intify = lambda x: (x[0], int(x[1]), int(x[2]))
    pat = re.compile('(.*):(.*)\.\.(.*)\(')
    map_info = (pd.read_table('gene_map_table_fb_2015_03.tsv', header=None,
                             names='symbol fbgn recomb cyto loc'.split(),
                             skiprows=6, index_col=1)
                .dropna(subset=['loc']))
    map_info = map_info.join(pd.DataFrame([intify(pat.match( map_info.ix[gene, 'loc']).groups())
                                           for gene in map_info.index],
                                          index=map_info.index,
                                          columns=['chr', 'start', 'stop']))

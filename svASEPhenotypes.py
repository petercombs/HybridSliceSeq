import pandas as pd
from collections import defaultdict

if __name__ == "__main__":
    logist = pd.read_table('analysis/results/logist.tsv', index_col=0)
    peak = pd.read_table('analysis/results/peak.tsv', index_col=0)
    fbgns = pd.read_table('prereqs/gene_map_table_fb_2016_01.tsv',
                          index_col=1, skiprows=5).ix[:, 0]

    svASE_genes = logist.index.union(peak.index)
    svASE_gns = tuple(fbgns[gene]+'[' for gene in svASE_genes)

    phenotypes = defaultdict(set)
    poor_phenotypes = ('viable', 'lethal', 'partially lethal', 'some die')
    for line in open('prereqs/allele_phenotypic_data_fb_2016_02.tsv'):
        line = line.split('\t')
        if (line and line[0].startswith(svASE_gns) and not
            line[2].startswith(poor_phenotypes)):
            phenotype = (line[2]
                         .split('with')[0]
                         .split('(')[0]
                         .split(',')[0]
                         .split('|')[0]
                         .strip()
                        )
            phenotypes[line[0].split('[')[0]].add(phenotype)






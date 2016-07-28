from __future__ import print_function
import pandas as pd
from collections import defaultdict

if __name__ == "__main__":
    logist = pd.read_table('analysis/results/logist.tsv', index_col=0)
    peak = pd.read_table('analysis/results/peak.tsv', index_col=0)
    fbgns = pd.read_table('prereqs/gene_map_table_fb_2016_01.tsv',
                          index_col=1, skiprows=5).ix[:, 0]

    svASE_genes = logist.index.union(peak.index)
    svASE_gns = tuple(fbgns[gene]+'[' for gene in svASE_genes)

    phenotypes_by_gene = defaultdict(set)
    genes_by_phenotype = defaultdict(set)
    poor_phenotypes = ('viable', 'lethal', 'partially lethal', 'some die')
    for line in open('prereqs/allele_phenotypic_data_fb_2016_02.tsv'):
        line = line.split('\t')
        if line and line[0].startswith(svASE_gns):
            if 'with' in line[2]:
                print(line[2])
                continue
            phenotype = (line[2]
                         #.split('with')[0]
                         .split('(')[0]
                         .split(',')[0]
                         .split('|')[0]
                         .strip()
                        )
            gene = line[0].split('[')[0]
            if not line[2].startswith(poor_phenotypes):
                phenotypes_by_gene[gene].add(phenotype)
            genes_by_phenotype[phenotype].add(gene)

    phenotypes_by_geneset = defaultdict(set)
    for phenotype in genes_by_phenotype:
        phenotypes_by_geneset[frozenset(genes_by_phenotype[phenotype])].add(phenotype)


    outfile = open('godot/results/phenotypes.dot', 'w')
    print("graph {", file=outfile)
    for gene in phenotypes_by_gene:
        for gene2 in phenotypes_by_gene:
            if gene == gene2: break
            intersection = phenotypes_by_gene[gene].intersection(phenotypes_by_gene[gene2])
            intersection.discard('viable')
            intersection.discard('visible')
            intersection.discard('fertile')
            if len(intersection) > 0:
                print('\t"{}" -- "{}" [weight={}, tooltip="{}"]; '
                      .format(gene, gene2, len(intersection),
                              ','.join(intersection)),
                      file=outfile)
    print("}", file=outfile)
    outfile.close()

    outfile = open('godot/results/pheno_genes.dot', 'w')

    print("graph {", file=outfile)
    print('\trankdir="LR"', file=outfile)
    for phenotype, genes in genes_by_phenotype.items():
        if (phenotype.startswith(('some die', 'long lived', 'short lived',
                                  'embryonic', 'larval', 'abdominal', 'embryo',
                                  'wild-type', 'male sterile', 'denticle',
                                  'neuroblast',
                                  'sterile', 'neuron', 'fertile', 'lethal', 'viable',
                                  'partially lethal', 'visible', ))
            or phenotype.endswith(('neuron', 'defective', 'segment'))):
            continue
        else:
            for gene in genes:
                print('\t"{}" -- "{}"; '
                      .format(phenotype, gene),
                     file=outfile)

    print("}", file=outfile)
    outfile.close()

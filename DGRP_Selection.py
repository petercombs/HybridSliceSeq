import pandas as pd

# Coordinates in R5 coordinates, 
# converted from FB listing and run through converter at 
# http://flybase.org/static_pages/downloads/COORD_625.html
target_genes = {
        'eloF':         ('3R', 5639883, 5640880),
        'trevelyan':    ('3R', 5638469, 5639661),
        'bond':         ('3R', 18376402, 18382216),
        'desatF':       ('3L', 11016635, 11017947),
        }



if 'dgrp' not in locals() or 'CHROM' not in locals()['dgrp'].columns:
    dgrp = (pd.read_table('Reference/dgrp_3.vcf', 
                          keep_default_na=False,
                          na_values=['./.'],
                          header=6,
                          )
           .dropna(how='all', axis=1)
           .rename(columns={"#CHROM":"CHROM"})
           )

target_SNPs = {}
target_SNPcounts = {}
for gene, (chrom, low, hi) in target_genes.items():
    tab =  dgrp.ix[(dgrp.CHROM == chrom) & (dgrp.POS > low-1000) & (dgrp.POS <  hi+1000)]
    target_SNPs[gene] = tab
    target_SNPcounts[gene] = pd.DataFrame([Counter(tab.ix[gene, 9:]) for gene in tab.index], index=tab.index)
    target_SNPcounts[gene]['CHR'] = tab.CHROM
    target_SNPcounts[gene]['POS'] = tab.POS
    target_SNPcounts[gene]['REF'] = tab.REF
    target_SNPcounts[gene]['ALT'] = tab.ALT


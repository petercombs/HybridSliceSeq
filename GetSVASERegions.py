import pandas as pd
from collections import Counter


if __name__ == "__main__":
    best_r2 = pd.read_table('analysis/results/svase_best', index_col=0,
                            header=None, squeeze=True)
    expr_nosvase = set(line.strip() for line in
                       open('analysis/results/expressed_nosvase.txt'))

    out_svase_bed = open('analysis/results/has_svase.bed', 'w')
    out_no_svase_bed = open('analysis/results/no_svase.bed', 'w')

    counts = Counter()
    for line in open('Reference/mel_r5_good.gtf'):
        line = line.split('\t')
        annot = {entr.split()[0]: entr.split()[1].strip('"')
                     for entr in line[-1].strip().strip(';').split('; ')}

        gene_name = annot['gene_name']
        bedline = "{}\t{}\t{}\t{}\n".format(
            line[0].strip('dmel_'), line[3], line[4], gene_name
        )
        if gene_name in best_r2 and best_r2[gene_name] > .5:
            out_svase_bed.write(bedline)
            counts['svase'] += 1
        elif gene_name in expr_nosvase:
            out_no_svase_bed.write(bedline)
            counts['nosvase'] += 1

    out_svase_bed.close()
    out_no_svase_bed.close()



import pandas as pd
from matplotlib.pyplot import (figure, violinplot, boxplot, xticks, savefig,
                               title, ylim, scatter)
from numpy.random import randn

if __name__ == "__main__":
    snpcount_kwargs = dict(
        header=None, names=['chr', 'start', 'stop', 'num_snps']
    )
    ratios_svase = {}
    ratios_nosvase = {}
    for featuretype in ['reg', 'tfbs', 'dnase', 'zld']:
        svASE_snpcounts = pd.read_table('analysis/results/has_svase_{}_snps.bed'
                                        .format(featuretype),
                                        **snpcount_kwargs)
        nosvASE_snpcounts = pd.read_table('analysis/results/no_svase_{}_snps.bed'
                                          .format(featuretype),
                                          **snpcount_kwargs)

        svASE_snpcounts['length'] = svASE_snpcounts.stop - svASE_snpcounts.start
        nosvASE_snpcounts['length'] = nosvASE_snpcounts.stop - nosvASE_snpcounts.start

        svASE_snpcounts['rate'] = svASE_snpcounts.num_snps / svASE_snpcounts.length
        nosvASE_snpcounts['rate'] = nosvASE_snpcounts.num_snps / nosvASE_snpcounts.length

        ratios_svase[featuretype] = svASE_snpcounts.rate
        ratios_nosvase[featuretype] = nosvASE_snpcounts.rate

        figure()
        violinplot([svASE_snpcounts.rate, nosvASE_snpcounts.rate], showmedians=True,
                   showextrema=False)
        boxplot([svASE_snpcounts.rate, nosvASE_snpcounts.rate])
        if featuretype != 'dnase':
            scatter(1 + .05 * randn(len(svASE_snpcounts)),
                    svASE_snpcounts.rate)
            scatter(2 + .05 * randn(len(nosvASE_snpcounts)),
                    nosvASE_snpcounts.rate)
        xticks([1, 2], ['svASE', 'no svASE'])
        ymin, ymax = ylim()
        ylim(-.01, ymax)
        title(featuretype)
        savefig('analysis/results/snp_rate_{}.png'.format(featuretype))



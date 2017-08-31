import pandas as pd
from matplotlib.pyplot import (figure, violinplot, boxplot, xticks, savefig,
                               title, ylim, scatter)
from numpy.random import randn
from os.path import exists

def get_datafiles(featuretype):
    svASE_snpcounts = pd.read_table('analysis/results/has_svase_{}_snpcounts_merged.bed'
                                    .format(featuretype),
                                    **snpcount_kwargs)
    if exists('analysis/results/no_svase_{}_snpcounts_merged.bed'.format(featuretype)):
        fname = 'analysis/results/no_svase_{}_snpcounts_merged.bed'.format(featuretype)
    elif exists('analysis/results/no_svase_{}_snpcounts.bed'.format(featuretype)):
        fname = 'analysis/results/no_svase_{}_snpcounts.bed'.format(featuretype)
    elif exists('analysis/results/no_svase_{}_clean_snpcounts.bed'.format(featuretype)):
        fname = 'analysis/results/no_svase_{}_clean_snpcounts.bed'.format(featuretype)
    else:
        raise ValueError("Can't find no svASE datafile")


    nosvASE_snpcounts = pd.read_table(fname,
                                          **snpcount_kwargs)

    svASE_snpcounts['length'] = abs(svASE_snpcounts.stop -
                                    svASE_snpcounts.start)
    nosvASE_snpcounts['length'] = abs(nosvASE_snpcounts.stop -
                                      nosvASE_snpcounts.start)

    svASE_snpcounts['rate'] = svASE_snpcounts.num_snps / svASE_snpcounts.length
    nosvASE_snpcounts['rate'] = nosvASE_snpcounts.num_snps / nosvASE_snpcounts.length
    return svASE_snpcounts, nosvASE_snpcounts


if __name__ == "__main__":
    snpcount_kwargs = dict(
        header=None, names=['chr', 'start', 'stop', 'num_snps']
    )
    ratios_svase = {}
    ratios_nosvase = {}
    svase_datasets = {}
    nosvase_datasets = {}
    for featuretype in ['reg', 'tfbs', 'dnase', 'zld', 'kr', 'med', 'D', 'da',
                        'twi' ]:
        svASE_snpcounts, nosvASE_snpcounts = get_datafiles(featuretype)
        svase_datasets[featuretype] = svASE_snpcounts
        nosvase_datasets[featuretype] = nosvASE_snpcounts


        ratios_svase[featuretype] = svASE_snpcounts.rate
        ratios_nosvase[featuretype] = nosvASE_snpcounts.rate

        figure()
        violinplot([svASE_snpcounts.num_snps.dropna(),
                    nosvASE_snpcounts.num_snps.dropna()], showmedians=True,
                   showextrema=False)
        boxplot([svASE_snpcounts.rate, nosvASE_snpcounts.rate])
        '''
        if featuretype != 'dnase':
            scatter(1 + .05 * randn(len(svASE_snpcounts)),
                    svASE_snpcounts.rate)
            scatter(2 + .05 * randn(len(nosvASE_snpcounts)),
                    nosvASE_snpcounts.rate)
                    '''
        xticks([1, 2], ['svASE', 'no svASE'])
        ymin, ymax = ylim()
        ylim(-.05 * ymax, ymax)
        title(featuretype)
        savefig('analysis/results/snp_rate_{}.png'.format(featuretype))



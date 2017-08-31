from scipy.stats import ttest_1samp
from scipy.stats import t as t_distribution
from numpy import array, isnan, nan, average, sqrt, inf
from progressbar import ProgressBar as pb
from collections import namedtuple
from argparse import ArgumentParser
import pandas as pd
import Utils as ut

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--min-genes', default=20, type=int)
    parser.add_argument('--test-only', default=False,
                        help="Test only go terms with a given code in second column")
    parser.add_argument('--drop-genes', default=None)
    parser.add_argument('--print-header', default=False, action='store_true')
    parser.add_argument('--print-counts', default=False, action='store_true')
    parser.add_argument('--pseudocounts', default=1, type=int)
    parser.add_argument('--smart-drop', default=False, action='store_true',
                        help='Use what we know about this dataset to remove '
                        'genes on the X in males')
    parser.add_argument('ase')
    parser.add_argument('categories')
    parser.add_argument('outfile')
    args = parser.parse_args()
    return args

def weighted_var(x, w, dropna=False):
    '''Translation of code from R to do wighted variance

    https://www.r-bloggers.com/weighted-t-test-in-r/

    # weighted variance, inspired by a function from Gavin Simpson on R-Help
    var.wt <- function(x, w, na.rm = FALSE) {
    if (na.rm) {
    w <- w[i <- !is.na(x)]
    x <- x[i]
    }
    sum.w <- sum(w)
    return((sum(w*x^2) * sum.w - sum(w*x)^2) / (sum.w^2 - sum(w^2)))
    }
    '''
    x = array(x)
    w = array(w)
    if dropna:
        i = ~(isnan(x) | isnan(w))
        x = x[i]
        w = w[i]
    ws = sum(w)
    return ( sum(w * (x**2)) * ws - sum(w * x)**2) / (ws**2 - sum(w**2))

weighted_ttest_result = namedtuple('weighted_ttest_result',
                                   ['estimate', 'se', 'conf_interval', 'tstat', 'df',
                                    'pvalue'])

def weighted_ttest(x, w, mu, conflevel=0.95, alternative='twosided', dropna=True):
    ''' Weighted t-test.

alternative: 'less', 'greater', 'twosided'

 weighted.t.test <- function(x, w, mu, conf.level = 0.95,
 alternative="two.sided", na.rm=TRUE) {

  if(!missing(conf.level) &
  (length(conf.level) != 1 || !is.finite(conf.level) ||
  conf.level < 0 || conf.level > 1))
  stop("'conf.level' must be a single number between 0 and 1")

   if (na.rm) {
   w <- w[i <- !is.na(x)]
   x <- x[i]
   }

    # to achieve consistent behavior in loops, return NA-structure in case of
    # complete missings
    if (sum(is.na(x)) == length(x))
        return(list(estimate=NA, se=NA, conf.int=NA, statistic=NA, df=NA, p.value=NA))

     # if only one value is present: this is the best estimate, no significance
     # test provided
     if (sum(!is.na(x)) == 1) {
     warning("Warning weighted.t.test: only one value provided; this value is
     returned without test of significance!", call.=FALSE)
     return(list(estimate=x[which(!is.na(x))], se=NA, conf.int=NA, statistic=NA,
 df=NA, p.value=NA))
 }

  x.w <- weighted.mean(x,w, na.rm=na.rm)
  var.w <- var.wt(x,w, na.rm=na.rm)
  df <- length(x)-1
  t.value <- sqrt(length(x))*((x.w-mu)/sqrt(var.w))
  se <- sqrt(var.w)/sqrt(length(x))

   if (alternative == "less") {
   pval <- pt(t.value, df)
   cint <- c(-Inf, x.w + se*qt(conf.level, df) )
   }
   else if (alternative == "greater") {
   pval <- pt(t.value, df, lower.tail = FALSE)
   cint <- c(x.w - se * qt(conf.level, df), Inf)
   }
   else {
   pval <- 2 * pt(-abs(t.value), df)
   alpha <- 1 - conf.level
   cint <- x.w + se*qt(1 - alpha/2, df)*c(-1,1)
   }

    names(t.value) <- "t"
    return(list(estimate=x.w, se=se, conf.int=cint, statistic=t.value, df=df,
p.value=pval))
    }
    '''
    x = array(x)
    w = array(w)
    if dropna:
        i = ~(isnan(x) | isnan(w))
        x = x[i]
        w = w[i]
    good_values = sum(~isnan(x))
    if good_values == 0:
        return weighted_ttest_result(estimate=nan, se=nan, conf_interval=nan,
                                     tstat=nan, df=nan, pvalue=nan)
    elif good_values == 1:
        # if only one value is present: this is the best estimate, no
        # significance test provided
        return weighted_ttest_result(estimate=x[~isnan(x)], se=nan, conf_interval=nan,
                          tstat=nan, df=nan, pvalue=nan)
    x_w = average(x, weights=w)
    var_w = weighted_var(x, w, dropna=dropna)
    df = len(x) - 1
    t_stat = sqrt(len(x)) * (x_w - mu) / sqrt(var_w)
    se = sqrt(var_w) / sqrt(len(x))
    if alternative == 'less':
        pval = t_distribution.cdf(t_stat, df)
        cint = (-inf, x_w + se * t_distribution.ppf(conflevel, df))
    elif alternative == 'greater':
        pval = t_distribution.sf(t_stat, df)
        cint = (x_w - se * t_distribution.ppf(conflevel, df), inf)
    else:
        pval = 2 * t_distribution.sf(abs(t_stat), df)
        alpha = 1-conflevel
        qt = t_distribution.ppf(1 - alpha/2, df)
        cint = (x_w - se * qt, x_w + se * qt)
    return weighted_ttest_result(estimate=x_w, se=se, conf_interval=cint, tstat=t_stat, df
                      = df, pvalue=pval)



if __name__ == "__main__":
    args = parse_args()
    drop_genes = set()
    if args.drop_genes:
        for gene in open(args.drop_genes):
            drop_genes.add(gene.strip())
    ase = (pd
           .read_table(args.ase, **ut.pd_kwargs)
           .select(**ut.sel_startswith(('melXsim', 'simXmel')))
           .rename(columns=lambda x: x if x.endswith('_ase_value') else x + '_ase_value')
           .rename(columns=lambda x: x.replace('melXsim_cyc14C_rep3',
                                               'melXsim_cyc14C_rep0'))
           .sort_index(axis=1)
          )
    ase.drop(drop_genes, inplace=True, errors='ignore')
    if args.smart_drop:
        chrom_of = ut.get_chroms()
        males = ('melXsim_cyc14C_rep3', 'simXmel_cyc14C_rep2')
        is_male  = [col.startswith(males) for col in ase.columns]
        ase.ix[chrom_of[ase.index] == 'X', is_male] = nan


    categories = {}
    for line in open(args.categories):
        line = line.strip().split()
        if ((len(line[2:]) < args.min_genes)
            or (args.test_only and line[1] != args.test_only)):
            continue
        categories[line[0]] = set(line[2:])

    tstats = pd.DataFrame(index=categories, columns=ase.columns, data=nan)
    pvals = pd.DataFrame(index=categories, columns=ase.columns, data=nan)
    ngenes = pd.DataFrame(index=categories, columns=ase.columns, data=nan)
    for sample in pb()(ase.columns):
        for category in categories:
            ase_vals = ase.ix[categories[category], sample].dropna()
            ngenes.ix[category, sample] = ase_vals.count()
            if ase_vals.count() < args.min_genes:
                continue
            t_result = ttest_1samp(ase_vals, 0)
            tstats.ix[category, sample] = t_result.statistic
            pvals.ix[category, sample] = t_result.pvalue

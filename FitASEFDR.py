from __future__ import print_function
import numpy as np
import pandas as pd
import sys
from os import getcwd
import os
from FitASEFuncs import (logistic, peak, fit_all_ase,
                         calculate_variance_explained)

from fyrd import Job
#import fyrd
#from multiprocessing import Pool
from progressbar import ProgressBar as pbar
from Utils import (sel_startswith, get_xs, pd_kwargs, get_chroms)
from queue import Queue
from argparse import ArgumentParser
from time import sleep

from warnings import filterwarnings
from pickle import dump

cluster_joblimit = 100
cluster_args = dict(time= '0:30:00', mem='60G',
                    partition='owners,hns,normal,hbfraser',
                    scriptpath='logs', outpath='logs', runpath=getcwd(),
                    cpus=4, cores=4)

def fit_and_eval(ase, func, xs, colnames, pool=False):
    res = fit_all_ase(ase, func, xs, colnames, pool, progress=True).dropna()
    return calculate_variance_explained( ase, xs, func, res)

def activate_job(waiting, active):
    if waiting.empty():
        return
    next_job = waiting.get()
    next_job.submit()
    active.put(next_job)


def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--expression', default='analysis_godot/summary.tsv',
                        type=str)
    parser.add_argument('--prefix', default='', type=str)
    parser.add_argument('--suffix', default='', type=str)
    parser.add_argument('data_to_fit',
                        default='analysis_godot/ase_summary_by_read.tsv',
                        type=str)
    return parser.parse_args()

if __name__ == "__main__":
    filterwarnings("ignore", ".*Covariance of the parameters.*",)
    filterwarnings("ignore", ".*overflow encountered in exp.*",)
    #expr = pd.read_table('analysis_godot/summary_fb.tsv', **pd_kwargs).dropna(how='all', axis=1)
    args = parse_args()
    ase = (pd
           .read_table(args.data_to_fit,
                       **pd_kwargs
                       )
           .dropna(how='all', axis=1)
           .dropna(how='all', axis=0)
           .select(**sel_startswith(('melXsim', 'simXmel')))
          )
    chrom_of = get_chroms()

    males = ('melXsim_cyc14C_rep3', 'simXmel_cyc14C_rep2')
    on_x = [chrom_of[gene] == 'X' if gene in chrom_of else False for gene in ase.index]
    is_male = [col.startswith(males) for col in ase.columns]
    ase.ix[on_x, is_male] = np.nan
    ase = ase.loc[ase.T.count() > len(ase.columns) / 2.0]

    hours = len(ase) / 1e4 * 1.5 + 2
    cluster_args['time'] = '{}:{}:00'.format(int(hours), int((hours % 1)*60))
    print("Estimate {} per iteration".format(cluster_args['time']))
    #cluster_args['queue'] = fyrd.Queue(user='self',
                                       #qtype=fyrd.queue.get_cluster_environment())
    print(cluster_args)
    sys.stdout.flush()


    xs = get_xs(ase)
    colnames = ['Amp', 'width', 'center', 'y_offset']
    peak_r2s = []
    logist_r2s = []

    n_perms = 1000
    waiting_jobs = Queue()
    active_jobs = Queue()
    for func, r2s in [(logistic, logist_r2s), (peak, peak_r2s)]:
        print('-'*30)
        print(func.__name__)
        print('-'*30, file=sys.stderr)
        print(func.__name__, file=sys.stderr)
        print('Building {} Jobs'.format(n_perms))
        sys.stdout.flush()
        sys.stderr.flush()
        for i in range(n_perms):
            print(i, end=' ')
            sys.stdout.flush()
            new_xs = pd.Series(index=xs.index, data=np.random.permutation(xs))
            waiting_jobs.put(Job(fit_and_eval,
                                 args=(ase, logistic, new_xs, colnames),
                                 kwargs={'pool': cluster_args['cpus']},
                                 suffix='_{}_{:04}'.format(func.__name__, i),
                                 **cluster_args
                                ))
            if i < cluster_joblimit:
                activate_job(waiting_jobs, active_jobs)

        sleep(60)
        for i in pbar(max_value=n_perms)(range(n_perms)):
            r2s.extend(active_jobs.get().get())
            if not waiting_jobs.empty():
                activate_job(waiting_jobs, active_jobs)

            dump({'logist': logist_r2s, 'peak': peak_r2s,
                  'last':'_{}_{:04}'.format(func.__name__, i), },
                 open('analysis/results/{prefix}fdr_{suffix}.pkl'
                      .format(prefix=args.prefix, suffix=args.suffix),
                      'wb'))

        np.save('analysis/results/{prefix}fdr_{name}{suffix}.numpy'
                             .format(prefix=args.prefix, name=func.__name__,
                                     suffix=args.suffix),
                np.array(r2s))

        np.array(r2s).tofile('analysis/results/{prefix}fdr_{name}{suffix}.numpy'
                             .format(prefix=args.prefix, name=func.__name__,
                                     suffix=args.suffix))

    np.array(peak_r2s).tofile('analysis/results/{prefix}fd_peak{suffix}.numpy'
                              .format(prefix=args.prefix, suffix=args.suffix))
    np.array(logist_r2s).tofile('analysis/results/{prefix}fd_logist{suffix}.numpy'
                                .format(prefix=args.prefix, suffix=args.suffix))

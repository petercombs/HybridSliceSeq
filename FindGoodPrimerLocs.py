from __future__ import print_function
from Bio import SeqIO
from subprocess import Popen
import pandas as pd
from collections import Counter, defaultdict
import re
from os import remove

startswith = lambda y: lambda x: x.startswith(y)

if __name__ == "__main__":
    parent = re.compile('FBgn[0-9]*')
    chrom = re.compile('loc=[^:]*')

    melix_orthologs = pd.read_table('prereqs/gene_orthologs_fb_2015_03.tsv',
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

    sim_trans = {parent.findall(rec.description)[0]:rec 
            for rec in SeqIO.parse('prereqs/dsim-all-translation-r2.01.fasta', 'fasta')
            if next(chrom.finditer(rec.description)).group().startswith('loc=Scf_X')}
    mel_trans = {rec.id:rec for rec in SeqIO.parse('prereqs/dmel-all-translation-r6.08.fasta', 'fasta')
            if next(chrom.finditer(rec.description)).group().startswith('loc=X')}

    fbgn_to_fbpp = defaultdict(list)
    for rec in mel_trans.values():
        fbgn = parent.findall(rec.description)[0]
        fbgn_to_fbpp[fbgn].append(rec.id)

    try:
        jobs = []
        for fbgn in fbgn_to_fbpp:
            if len(fbgn_to_fbpp[fbgn]) > 1: continue
            if fbgn not in melix_orthologs.index: continue
            sim_fbgn = melix_orthologs.ix[fbgn, 'Ortholog_FBgn_ID']
            if not isinstance(sim_fbgn, str): continue
            if sim_fbgn not in sim_trans: continue

            SeqIO.write(mel_trans[fbgn_to_fbpp[fbgn][0]], 
                    'tmp/{}.fasta'.format(fbgn),
                    'fasta')
            SeqIO.write(sim_trans[sim_fbgn],
                    'tmp/{}.fasta'.format(sim_fbgn),
                    'fasta')
            jobs.append(Popen(['blat', '-prot',
                '-noHead',
                'tmp/{}.fasta'.format(fbgn),
                'tmp/{}.fasta'.format(sim_fbgn),
                'tmp/{}-{}.psl'.format(fbgn, sim_fbgn),
                ],
                stdout=open('/dev/null', 'w'),
                stderr=open('/dev/null', 'w'),
                ))
        for job in jobs:
            finish = job.wait()
            if finish != 0: continue
            res = pd.read_table(job.args[-1],
                    #skiprows=5,
                    usecols=[5,7,9,13],
                    names=['Q_gaps', 'T_gaps', 'Q_name', 'T_name']
                    )
            if len(res) == 0: continue
            if max(abs(res.Q_gaps - res.T_gaps)) > 10:
                print('\t'.join(str(i) for i in res.ix[0]))
            else:
                try:
                    remove(job.args[-1])
                    remove(job.args[-2])
                    remove(job.args[-3])
                except OSError:
                    print("!", end='')

    finally:
        for job in jobs:
            job.wait()

            

            

        


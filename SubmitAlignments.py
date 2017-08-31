from subprocess import Popen, PIPE
from argparse import ArgumentParser
from progressbar import ProgressBar as pb
from io import StringIO
import itertools

parser = ArgumentParser()
parser.add_argument('species1', type=str, 
        help='The first species (e.g. "mel")')
parser.add_argument('species2', type=str,
        help='The second species (e.g. "sim")')

args = parser.parse_args()

script = open('qsub_base.sh').read()

species1_gtf = open('Reference/{}_good.gtf'.format(args.species1))
species2_gtf = open('Reference/{}_good.gtf'.format(args.species2))

species1_chroms = {line.split()[0] for line in species1_gtf}
species2_chroms = {line.split()[0] for line in species2_gtf}

print(species1_chroms, species2_chroms)
items = list(itertools.product(species1_chroms, species2_chroms))

jobs = []
for id1, id2 in pb()(items):
    job = script.format(
            job_name=id1+'_'+id2,
            id1 = id1,
            id2 = id2, 
            species1=args.species1,
            species2=args.species2)
    jobs.append(Popen(['qsub'], stdin=PIPE))
    jobs[-1].communicate(bytes(job, 'ASCII'))

for job in jobs:
    job.wait()


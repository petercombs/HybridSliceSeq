from subprocess import Popen
from tempfile import NamedTemporaryFile

script = open('qsub_base.sh').read()

species1_gtf = open('Reference/mel_good.gtf')
species2_gtf = open('Reference/sim_good.gtf')

species1_chroms = {line.split()[0] for line in species1_gtf}
species2_chroms = {line.split()[0] for line in species2_gtf}

print(species1_chroms, species2_chroms)

for id1 in species1_chroms:
    for id2 in species2_chroms:
        outfile = NamedTemporaryFile(dir='.')
        job = script.format(
            job_name=id1+'_'+id2,
            id1 = id1,
            id2 = id2)
        print(job)
        outfile.file.write(bytes(job, 'UTF-8'))
        outfile.flush()
        Popen(['qsub', outfile.name]).wait()
        outfile.close()
        


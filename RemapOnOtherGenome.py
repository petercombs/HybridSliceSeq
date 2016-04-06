from pysam import Samfile
from sys import argv
from subprocess import Popen

if __name__ == "__main__":
    infile = Samfile(argv[1])
    on_target = argv[2]
    cl = infile.header['PG'][0]['CL'].split()
    i = cl.index('--genomeDir')
    cl[i+1] = on_target.replace('/Genome', '')
    print(' '.join(cl))
    out = Popen(cl)
    out.wait()
    if cl[0] == 'STAR':
        samfile = cl[cl.index('--outFileNamePrefix')+1]+'/Aligned.out.sam'
        Popen(['samtools', 'view', '-bS', '-o', argv[3],
            samfile,
            ]).wait()
        Popen(['rm', samfile]).wait()





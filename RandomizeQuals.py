""" 
Randomize Quality Scores

When using samtools rmdup, only the highest quality read is kept, which is not
what we want to do for ASE.  Therefore, one solution is to randomize the
quality score for each read, then feed it through to samtools, then later
restore the correct quality score (or not).  It takes a few more passes through
the data than WASP, but unlike WASP, it's simple and it actually works.  

Sample Usage: 

    python RandomizeQuals.py file.bam - | samtools rmdup - - | python RestoreQuals.py - out.bam

"""
from __future__ import print_function
from pysam import Samfile
from random import randint
from sys import argv, exit

if __name__ == "__main__":
    if len(argv) != 3:
        print("""
Randomize Map Quality Scores

Usage:
python {} INFILE OUTFILE
    INFILE: SAM or BAM formatted file (or - for stdin)
    OUTFILE: File to write to (SAM or BAM or - for BAM-formmated stdout)
    """.format(argv[0]))
        exit(1)
    with Samfile(argv[1]) as insam:
        mode = 'w' if argv[2].endswith('sam') else 'wb'
        outsam = Samfile(argv[2], mode, template=insam)
        for read in insam:
            read.set_tag('oq', read.mapq)
            read.mapq = randint(0,254)
            outsam.write(read)

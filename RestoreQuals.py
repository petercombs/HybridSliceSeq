from pysam import Samfile
from sys import argv, exit
"""
Restore Quality Scores

"""

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
            read.mapq = read.get_tag('oq')
            outsam.write(read)

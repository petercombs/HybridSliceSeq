from pysam import Samfile
from os import path
from sys import argv

sf = Samfile(argv[1])

of = open(path.join(path.dirname(argv[1]), path.basename(argv[1]) + '.chromsizes'), 'w')
print(of.name)
for rec in sf.header['SQ']:
    if 'SN' in rec:
        of.write('{}\t{}\n'.format(rec['SN'], rec['LN']))


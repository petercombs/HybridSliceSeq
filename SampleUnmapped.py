from sys import argv
from random import randint

max_els = int(1e4)

if __name__ == "__main__":
    out = []
    infile = open(argv[1])
    k = 0
    for k, line in enumerate(infile):
        if k % 4 == 1:
            out.append(line.strip())
        if len(out) >= max_els:
            break
    
    next(infile)
    next(infile)
    k //= 4
    for i, line in enumerate(infile):
        if i % 4 == 1:
            j = randint(0, i//4 + k)
            if j < max_els:
                out[j] = line.strip()

    with open(argv[2], 'w') as outf:
        for i, el in enumerate(out):
            outf.write('>sample_{:06}\n{}\n'.format(i+1, el))





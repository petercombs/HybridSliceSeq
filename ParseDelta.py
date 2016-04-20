from numpy import (zeros, arange, sum, array, argwhere, abs)

def parse_records(delta_file):
    while True:
        try:
            yield parse_record(delta_file)
        except StopIteration:
            break



def parse_record(delta_file, return_names=True):
    line = next(delta_file)
    while not line.startswith('>'):
        line = next(delta_file)
    # The first line is
    #  >seq1 seq2 length1 length2
    line = line.strip().strip('>').split()
    n1 = line[0]
    n2 = line[1]
    l1 = int(line[-2])
    l2 = int(line[-1])
    # The next line is
    # start1 end1 start2 end2 ??? ???
    line = next(delta_file).strip().split()
    i1 = int(line[0])-1
    p1 = i1
    e1 = int(line[1])
    i2 = int(line[2])-1
    e2 = int(line[3])
    p2 = i2

    offset = i2 - i1
    deltas = []
    while line:
        line = next(delta_file)
        line = int(line)
        deltas.append(line)
    deltas = array(deltas)

    longest = max(l1 + sum(deltas < 0) - offset, l2 + sum(deltas > 0)+offset)
    pos1 = zeros(longest, dtype=int)
    pos2 = zeros(longest, dtype=int)
    offset1 = offset if offset > 0 else 0
    offset2 = -offset if offset < 0 else 0
    i1 += offset1
    i2 += offset2
    pos1[offset1:i1] = arange(p1)
    pos2[offset2:i2] = arange(p2)

    for d in deltas:
        if d == 0:
            break
        ad = abs(d)
        pos1[i1:i1+ad] = arange(p1, p1+ad)
        pos2[i2:i2+ad] = arange(p2, p2+ad)
        if d > 0:
            p1 += ad
            p2 += ad
            i1 += ad
            i2 += ad
            p2 -= 1
        elif d < 0:
            p1 += ad
            p2 += ad
            i1 += ad
            i2 += ad
            p1 -= 1
        else:
            assert False




    pos1[i1:i1+l1-p1] = arange(p1, l1)
    pos2[i2:i2+l2-p2] = arange(p2, l2)
    pos1[i1+l1-p1:] = pos1[i1+l1-p1-1]
    pos2[i2+l2-p2:] = pos2[i2+l2-p2-1]


    if return_names:
        return (n1, pos1), (n2, pos2)
    return pos1, pos2


if __name__ == "__main__":
    pass




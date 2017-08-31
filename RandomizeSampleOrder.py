from __future__ import print_function
import random
import itertools

random.seed(0)

#samples_MA = ['MA{:02}'.format(n) for n in range(35,63)]

samples = ['{}{:02}'.format(group, n)
           for group, first, last in [
               ('MB', 2, 30),     # 50%
               ('MC', 4, 34),     # 47%
               #('MD', 3, 26),     # 44%
               ('SB', 6, 34),     # 51%
               #('RB', 5, 31),     # 58%
               #('RC', 14, 19),    # 59%
               ('RD', 7, 33),     # 55%
               ('RE', 3, 31),     # 52%
           ]
           for n in range(first, last+1)
           ]

random.shuffle(samples)

samps = (s for s in samples)
print('Strip\tTube\tSample')
strip_number = 1
for tube_numbers in itertools.repeat(range(1, 9)):
    for tube in tube_numbers:
        s = next(samps)
        print(strip_number, tube, s, sep='\t')
    strip_number += 1


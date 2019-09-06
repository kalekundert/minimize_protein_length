#!/usr/bin/env python3

import numpy as np
from pymol import cmd

scores = np.fromfile('deletion_scores.dat')

low = np.percentile(scores, 10)
high = np.percentile(scores, 90)

sele = 'chain B and polymer.protein'
cmd.alter(sele, 'b=scores[int(resi)]')
cmd.spectrum('b', 'blue_white_red', sele, minimum=low, maximum=high)

print(f'low:  {low:.2f}')
print(f'high: {high:.2f}')





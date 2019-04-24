#!/usr/bin/env python

import sys, os, re
import numpy as np

np.set_printoptions(precision=3, suppress=True, threshold=np.inf, linewidth=300)

TMP_DIR = 'tmp'
INDEX = 'v2015/INDEX.csv'

onames = {
    'X': 'matlab_X.txt',
    'y': 'matlab_y.txt',
    'code': 'matlab_code.txt',
}

value_map = {}
for i, line in enumerate(open(INDEX, 'rt')):
    if i == 0:
        continue
    it = line.rstrip().split(',')
    code = it[0]
    val = float(it[1])
    value_map[code] = val


pat = re.compile('/([1-9][0-9a-z][0-9a-z][0-9a-z])/')
Xs = []
codes = []
ys = []
count = 0
for name in os.listdir(TMP_DIR):
    if name.endswith('_out.npz'):
        fname = f'{TMP_DIR}/{name}'
        d = np.load(fname)
        pro_fname = d['pro_fname'].item()
        m = pat.search(pro_fname)
        if m:
            code = m.group(1)
        else:
            continue
        y = value_map.get(code)
        if y is None:
            continue
        count += 1
        print(count, name, file=sys.stderr)
        P = d['P']
        Xs.append(P)
        ys.append(y)
        codes.append(code)

Xs = np.array(Xs)
ys = np.array(ys)
codes = np.array(codes)
print(Xs.shape, file=sys.stderr)

np.savez('Xy.npz', X=Xs, codes=codes, y=ys)

### FOR MATLAB

matlab_X = Xs.reshape((Xs.shape[0], Xs.shape[1]*Xs.shape[2]*Xs.shape[3])).T
print(matlab_X.shape, file=sys.stderr)

out = open(onames['X'], 'wt')
for i in range(matlab_X.shape[0]):
    print(' '.join(['%.3f' % v for v in matlab_X[i].tolist()]), file=out)
out.close()

out = open(onames['y'], 'wt')
print(' '.join(['%.3f' % v for v in ys.tolist()]), file=out)
out.close()

out = open(onames['code'], 'wt')
print(' '.join(codes.tolist()), file=out)
out.close()

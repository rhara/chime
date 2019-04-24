#!/usr/bin/env python

import sys, os, shlex
import multiprocessing as mp
import subprocess as sp

here = os.path.dirname(os.path.abspath(__file__))

srcdir = sys.argv[1]

LIGAND_ATOMS = 90

def worker(name):
    lig_fname = '%s/%s/%s_ligand.sdf' % (srcdir, name, name)
    pro_fname = '%s/%s/%s_pocket.pdb.gz' % (srcdir, name, name)
    if os.path.exists(lig_fname) and \
       os.path.exists(pro_fname):
        cmd = f'python {here}/atomic_conv.py {lig_fname} {pro_fname} tmp/{name}_out.npz --padding {LIGAND_ATOMS}'
        proc = sp.Popen(shlex.split(cmd), universal_newlines=True, stdout=sp.PIPE, stderr=sp.STDOUT)
        # out = proc.stdout.read().rstrip()
        # print(out)
        proc.wait()
        ret = True
    else:
        ret = False
    return name, ret

os.makedirs('tmp', exist_ok=True)
for name in os.listdir('tmp'):
    os.unlink(f'tmp/{name}')

ls = sorted([name for name in os.listdir(srcdir) if len(name) == 4 and name[0] in '0123456789'])
pool = mp.Pool(mp.cpu_count())
try:
    count = 0
    for name, ret in pool.imap_unordered(worker, ls):
        if ret:
            count += 1
            print(count, name)
except KeyboardInterrupt:
    sys.exit(1)

#!/usr/bin/env python

import sys
import numpy as np
from chime import Timer, ComplexFeaturizer, readLigand, readProtein

np.set_printoptions(precision=3, linewidth=300)


class Opts:
    """
    Option parser specific to this script
    """
    def __init__(self, argv):
        argv = list(argv)
        self.neighbor_size = 12
        self.radials_setup = '1.5 6.0  0.5'
        self.padding = 90
        args = []
        while argv:
            a = argv.pop(0)
            if a in ['--neighbor-size', '-M']:
                self.neighbor_size = int(argv.pop(0))
                continue
            if a in ['--radials-setup', '-R']:
                self.radials_setup = argv.pop(0)
                continue
            if a in ['--padding', '-P']:
                self.padding = int(argv.pop(0))
                continue
            args.append(a)
        self.prog_name = args[0]
        self.lig_fname = args[1]
        self.pro_fname = args[2]
        self.out_fname = args[3]


def main(argv):
    opts = Opts(argv)
    radials = tuple(map(float, opts.radials_setup.split()))

    lig = readLigand(opts.lig_fname, sanitize=False)
    lig.SetProp('fname', opts.lig_fname)
    pro = readProtein(opts.pro_fname)
    pro.SetProp('fname', opts.pro_fname)

    if not lig or opts.padding < lig.GetNumAtoms():
        natoms = 0
        try:
            natoms = lig.GetNumAtoms()
        except:
            pass
        print(f'Error "{opts.lig_fname}" natoms={natoms} is too big and is skipped', file=sys.stderr)
        sys.exit(-1)

    featurizer = ComplexFeaturizer(atomtype_list=[6,7,8,16],
                                   M=opts.neighbor_size,
                                   radials=radials,
                                   padding=opts.padding)

    print('radials=%s' % featurizer.radials, file=sys.stderr)

    timer = Timer()
    featurizer(lig, pro)
    print('%.2fs' % timer.elapsed(), file=sys.stderr)

    R, Z, E, P = featurizer.R, featurizer.Z, featurizer.E, featurizer.P

    print('Shapes: R:%s, Z:%s, E:%s, P:%s' % (R.shape, Z.shape, E.shape, P.shape), file=sys.stderr)

    featurizer.save(opts.out_fname)


if __name__ == '__main__':
    main(sys.argv)

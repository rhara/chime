import sys
import numpy as np
from chime import Timer, ComplexFeaturizer

np.set_printoptions(precision=3, linewidth=300)


class Opts:
    """
    Option parser specific to this script
    """
    def __init__(self, argv):
        argv = list(argv)
        self.neighbor_size = 12
        self.radials_setup = '1.5 12.0  0.5'
        args = []
        while argv:
            a = argv.pop(0)
            if a in ['--neighbor-size', '-M']:
                self.neighbor_size = int(argv.pop(0))
                continue
            if a in ['--radials-setup', '-R']:
                self.radials_setup = argv.pop(0)
                continue
            args.append(a)
        self.prog_name = args[0]
        self.lig_fname = args[1]
        self.pro_fname = args[2]
        self.out_fname = args[3]


def main(argv):
    opts = Opts(argv)
    rmin, rmax, rint = map(float, opts.radials_setup.split())
    radials = (rmin, rmax, rint)

    featurizer = ComplexFeaturizer(M=opts.neighbor_size,
                                   radials=radials)

    print('radials=%s' % featurizer.radials, file=sys.stderr)

    timer = Timer()
    featurizer(opts.lig_fname, opts.pro_fname)
    print('%.2fs' % timer.elapsed(), file=sys.stderr)

    R, Z, E, P = featurizer.R, featurizer.Z, featurizer.E, featurizer.P

    print('Shapes: R:%s, Z:%s, E:%s, P:%s' % (R.shape, Z.shape, E.shape, P.shape), file=sys.stderr)

    featurizer.save(opts.out_fname)


if __name__ == '__main__':
    main(sys.argv)

import sys
from rdkit import Chem
# import numpy as np
from lmgcn.matrix import getDistanceMatrix, getAdjacentMatrix
from lmgcn.chemio import readLigand, readProtein
from lmgcn.timer import Timer

class Opts:
    """
    Option parser specific to this script
    """
    def __init__(self, argv):
        argv = list(argv)
        self.verbose = False
        self.args = []
        while argv:
            a = argv.pop(0)
            if a == '--verbose':
                self.verbose = True
                continue
            self.args.append(a)
        self.progname = self.args[0]
        self.ligname = self.args[1]
        self.protname = self.args[2]

def main(argv):
    opts = Opts(argv)

    ligand = readLigand(opts.ligname)
    protein = readProtein(opts.protname)
    complex_mol = Chem.CombineMols(ligand, protein)

    timer = Timer()
    R = getDistanceMatrix(complex_mol)
    A = getAdjacentMatrix(complex_mol, bondorder=False)
    print('%.2fs' % timer.elapsed(), file=sys.stderr)

    if opts.verbose:
        print(R)
        print(A)

if __name__ == '__main__':
    main(sys.argv)

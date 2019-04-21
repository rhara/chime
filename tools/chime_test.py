#!/usr/bin/env python

import sys
import numpy as np
from rdkit import Chem
from chime import (
    Timer,
    readLigand,
    readProtein,
    getAtomPairGraphDistanceMatrix,
    getActivePocket,
)

np.set_printoptions(precision=3, edgeitems=10, linewidth=300)


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
    print(ligand.GetNumAtoms())
    print(protein.GetNumAtoms())
    complex_mol = Chem.CombineMols(ligand, protein)

    timer = Timer()

    gdm = getAtomPairGraphDistanceMatrix(ligand)
    print(gdm)

    pocket = getActivePocket(ligand, protein)
    print('%.2fs' % timer.elapsed(), file=sys.stderr)

    print(pocket.GetNumAtoms())

    writer = Chem.PDBWriter('out.pdb')
    writer.write(pocket)
    writer.close()

if __name__ == '__main__':
    main(sys.argv)

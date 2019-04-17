import sys, gzip
from rdkit import Chem
import numpy as np

def getDistanceMatrix(mol):
    coords = mol.GetConformer(0).GetPositions()

    natoms = mol.GetNumAtoms()
    matrix = np.zeros((natoms, natoms))
    bonds = np.zeros((natoms, natoms), dtype=np.int)

    for i in range(natoms):
        cs = np.tile(coords[i], natoms).reshape((natoms, 3))
        d = np.linalg.norm(coords - cs, axis=1)
        matrix[i] = d

    return matrix

def getAdjacentMatrix(mol):
    natoms = mol.GetNumAtoms()
    matrix = np.zeros((natoms, natoms))
    for i in range(natoms-1):
        for j in range(i+1, natoms):
            bond = mol.GetBondBetweenAtoms(i, j)
            if bond:
                matrix[i, j] = matrix[j, i] = 1

    return matrix


class Opts:
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


opts = Opts(sys.argv)

ligand = Chem.MolFromMolBlock((gzip.open if opts.ligname.endswith('.gz') else open)(opts.ligname, 'rt').read())
protein = Chem.MolFromPDBBlock((gzip.open if opts.protname.endswith('.gz') else open)(opts.protname, 'rt').read())
complex_mol = Chem.CombineMols(ligand, protein)

R = getDistanceMatrix(complex_mol)
A = getAdjacentMatrix(complex_mol)

if opts.verbose:
    print(R)
    print(A)

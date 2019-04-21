"""
Reference

Molecular Graph Convolutions: Moving Beyond Fingerprints
S. Kearnes, K. McCloskey, M. Berndl, V. Pande, P. Riley
arXiv:1603.00856v3 (2016)

Atom features:
    Atom type
    Chirality
    Formal charge
    Parital charge
    Ring size
    Hybridization
    Hydrogen bonding
    Aromaticity

Atom pair features:
    Bond type
    Graph distance
    Same ring


"""

from rdkit import Chem
from .atom import getAtomRingSize
import numpy as np


__donor_pat = ['$([N;!H0;v3,v4&+1])',
               '$([O,S;H1;+0])',
               'n&H1&+0']

__acceptor_pat = ['$([O,S;H1;v2;!$(*-*=[O,N,P,S])])',
                  '$([O,S;H0;v2])',
                  '$([O,S;-])',
                  '$([N;v3;!$(N-*=[O,N,P,S])])',
                  'n&H0&+0',
                  '$([o,s;+0;!$([o,s]:n);!$([o,s]:c:n)])']

DONOR_SMARTS = Chem.MolFromSmarts('[' + ','.join(__donor_pat) + ']')
ACCEPTOR_SMARTS = Chem.MolFromSmarts('[' + ','.join(__acceptor_pat) + ']')


def getAtomTypeVector(mol):
    """
    Extension of Atom type
    """
    return np.array([a.GetAtomicNum() for a in mol.GetAtoms()], dtype=np.uint8)


def getAtomChiralityVector(mol):
    """
    Atom Chirality
    """
    return [int(a.GetChiralTag()) for a in mol.GetAtoms()]


def getAtomFormalChargeVector(mol):
    """
    Atom Formal charge
    """
    return np.array([a.GetFormalCharge() for a in mol.GetAtoms()], dtype=np.int8)


def getAtomPartialChargeVector(mol):
    """
    Atom Partial charge
    """
    from rdkit.Chem import rdPartialCharges
    rdPartialCharges.ComputeGasteigerCharges(mol)
    return np.array([float(a.GetProp('_GasteigerHCharge')) for a in mol.GetAtoms()], dtype=np.float16)


def getAtomRingSizeVector(mol):
    """
    Atom Ring sizes
    """
    return np.array([getAtomRingSize(a) for a in mol.GetAtoms()], dtype=np.int16)


def getAtomHybridizationVector(mol):
    """
    Atom Hybridization
    """
    return [int(a.GetHybridization()) for a in mol.GetAtoms()]


def getAtomHydrogenBondingVector(mol):
    """
    Atom Hydrogen bonding
    (note: not one-hot expression)
    donor 1
    acceptor 2
    """
    N = mol.GetNumAtoms()
    v = [0]*N
    donors = mol.GetSubstructMatches(DONOR_SMARTS)
    acceptors = mol.GetSubstructMatches(ACCEPTOR_SMARTS)
    for donor in donors:
        v[donor[0]] |= 1
    for acceptor in acceptors:
        v[acceptor[0]] |= 2

    return np.array(v)


def getAtomAromaticVector(mol):
    """
    Atom Aromaticity
    """
    return np.array([a.GetIsAromatic() for a in mol.GetAtoms()], dtype=np.uint8)


def getAtomPairBondTypeMatrix(mol, bondorder=False):
    """
    Extension of Atom pair Bond type matrix (symmatric)
    if bondorder == True, then bondtype matrix
    if bondorder == False, then adjacent matrix
    """
    natoms = mol.GetNumAtoms()
    matrix = np.zeros((natoms, natoms), dtype=np.float16 if bondorder else np.uint8)
    for i in range(natoms-1):
        for j in range(i+1, natoms):
            bond = mol.GetBondBetweenAtoms(i, j)
            if bond:
                if bondorder:
                    matrix[i, j] = matrix[j, i] = bond.GetBondTypeAsDouble()
                else:
                    matrix[i, j] = matrix[j, i] = 1

    return matrix


def getAtomPairGraphDistanceMatrix(mol, maxlen=0):
    """
    Atom pair Graph distance matrix
    (BFS)
    """
    def rec(idxs):
        if len(idxs) == 0:
            return
        depths.append(idxs)
        for idx in idxs:
            visited.add(idx)
        if maxlen and maxlen < len(depths):
            return
        nidxs = set()
        for idx in idxs:
            for nidx in np.where(A[idx] == 1)[0]:
                nidx = nidx.item()
                if nidx not in visited:
                    nidxs.add(nidx)
        rec(nidxs)

    N = mol.GetNumAtoms()
    A = getAtomPairBondTypeMatrix(mol)

    mat = np.zeros((N, N), dtype=np.uint16)
    depths = []
    visited = set()
    for i in range(N):
        depths.clear()
        visited.clear()
        rec(set([i]))
        for j in range(1, len(depths)):
            for k in depths[j]:
                mat[i, k] = mat[k, i] = j
    return mat


def getAtomPairSameRingMatrix(mol):
    """
    Atom pair Same ring matrix
    """

    def insamering(mol, i, j):
        paths = []
        def rec(i, goal, visited=[]):
            atom = mol.GetAtomWithIdx(i)
            if not atom.IsInRing():
                return
            if i == goal:
                paths.append(visited + [i])
                return
            if i in visited:
                return
            visited.append(i)
            for neighbor in atom.GetNeighbors():
                nidx = neighbor.GetIdx()
                rec(nidx, goal, list(visited))

        rec(i, j)
        return 1 if paths else 0

    N = mol.GetNumAtoms()
    ap = np.zeros((N, N), dtype=np.uint16)

    for i in range(N-1):
        atom1 = mol.GetAtomWithIdx(i)
        if not atom1.IsInRing():
            continue
        for j in range(i+1, N):
            atom2 = mol.GetAtomWithIdx(j)
            if not atom2.IsInRing():
                continue
            ap[i, j] = ap[j, i] = insamering(mol, i, j)

    return ap


def getAtomPairDistanceMatrix(mol1, mol2=None, dtype=np.float16):
    if mol2 is None:
        mol2 = mol1
    N1 = mol1.GetNumAtoms()
    N2 = mol2.GetNumAtoms()
    xyzs1 = mol1.GetConformer(0).GetPositions()
    xyzs2 = mol2.GetConformer(0).GetPositions()
    matrix = np.zeros((N1, N2), dtype=dtype)

    for i in range(N1):
        cs = np.tile(xyzs1[i], N2).reshape((N2, 3))
        matrix[i] = np.linalg.norm(xyzs2 - cs, axis=1)

    return matrix


def getNeighborList(mol1, mol2=None, M=6):
    dist_mat = getAtomPairDistanceMatrix(mol1, mol2)
    target = mol2 if mol2 else mol1
    atomtype_array = getAtomTypeVector(target)
    N = mol1.GetNumAtoms()
    nbors = np.zeros([N, M], dtype=np.int16)
    R = np.zeros([N, M], dtype=np.float16)
    Z = np.zeros([N, M], dtype=np.int16)
    for i in range(N):
        nbors[i] = np.argsort(dist_mat[i])[:M]
        R[i] = dist_mat[i, nbors[i]]
        Z[i] = atomtype_array[nbors[i]]

    return R, Z

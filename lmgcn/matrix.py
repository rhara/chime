from rdkit import Chem
import numpy as np

def getAtomicNumVector(mol):
    return np.array([a.GetAtomicNum() for a in mol.GetAtoms()], dtype=np.uint8)


def getDistanceMatrix(mol1, mol2=None, dtype=np.float16):
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

def getAdjacentMatrix(mol, bondorder=False):
    """
    Calculate adjacent matrix of np.array
    if bondorder, returns matrix of np.float16, else returns that of np.int

    Returns symmetric matrix
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

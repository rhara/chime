from rdkit import Chem
import numpy as np

def getDistanceMatrix(mol, dtype=np.float16):
    """
    Calculate distance matrix of np.array
    """
    coords = mol.GetConformer(0).GetPositions()

    natoms = mol.GetNumAtoms()
    matrix = np.zeros((natoms, natoms), dtype=dtype)

    for i in range(natoms):
        cs = np.tile(coords[i], natoms).reshape((natoms, 3))
        d = np.linalg.norm(coords - cs, axis=1)
        matrix[i] = d

    return matrix

def getAdjacentMatrix(mol, bondorder=False):
    """
    Calculate adjacent matrix of np.array
    if bondorder, returns matrix of np.float16, else returns that of np.int

    Returns symmetric matrix
    """
    natoms = mol.GetNumAtoms()
    matrix = np.zeros((natoms, natoms), dtype=np.float16 if bondorder else np.int8)
    for i in range(natoms-1):
        for j in range(i+1, natoms):
            bond = mol.GetBondBetweenAtoms(i, j)
            if bond:
                if bondorder:
                    matrix[i, j] = matrix[j, i] = bond.GetBondTypeAsDouble()
                else:
                    matrix[i, j] = matrix[j, i] = 1

    return matrix



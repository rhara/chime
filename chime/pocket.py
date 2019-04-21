from .graph_features import getAtomPairDistanceMatrix
from .submol import getSubmol
import numpy as np
from itertools import product


def getActivePocket(ligand, protein, distance=6.0, minclusize=5, gap=5):
    dmat = getAtomPairDistanceMatrix(ligand, protein) < distance

    vres = np.zeros(protein.GetNumAtoms(), dtype=np.uint16)
    for i in range(protein.GetNumAtoms()):
        vres[i] = protein.GetAtomWithIdx(i).GetPDBResidueInfo().GetResidueNumber()

    res_valid = set()
    for i, j in product(range(dmat.shape[0]), range(dmat.shape[1])):
        if dmat[i, j]:
            res_valid.add(vres[j])
    res_valid = sorted(res_valid)

    clusters = []
    res_tmp = set()
    for r in res_valid:
        if clusters and r - clusters[-1][-1] == 1:
            clusters[-1].append(r)
        else:
            clusters.append([r])
    for clu in clusters:
        while len(clu) < minclusize:
            clu.insert(0, clu[0]-1)
            clu.append(clu[-1]+1)
        for r in clu:
            res_tmp.add(r)
    res_valid = sorted(res_tmp)

    aidxs = []
    for i in range(protein.GetNumAtoms()):
        if vres[i] in res_valid:
            aidxs.append(i)
    aidxs = sorted(aidxs)

    return getSubmol(protein, aidxs, confs=True)

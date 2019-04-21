from rdkit import Chem
from itertools import combinations

def getSubmol(mol, aidxs, confs=False):
    submol = Chem.RWMol()
    for aix in aidxs:
        submol.AddAtom(mol.GetAtomWithIdx(aix))
    for i, j in combinations(aidxs, 2):
        bond = mol.GetBondBetweenAtoms(i, j)
        if bond:
            submol.AddBond(aidxs.index(i),
                           aidxs.index(j),
                           bond.GetBondType())
    if confs:
        for conf in mol.GetConformers():
            newconf = Chem.Conformer(len(aidxs))
            for i in range(len(aidxs)):
                newconf.SetAtomPosition(i, conf.GetAtomPosition(aidxs[i]))
                newconf.SetId(conf.GetId())
                newconf.Set3D(conf.Is3D())
            submol.AddConformer(newconf)
    return submol

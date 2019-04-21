import gzip
from rdkit import Chem


def readLigand(fname, **kwargs):
    mol = Chem.MolFromMolBlock((gzip.open if fname.endswith('.gz') else open)(fname, 'rt').read(), **kwargs)
    # Chem.AssignAtomChiralTagsFromStructure(mol)
    # Chem.AssignStereochemistry(mol)
    # Chem.SanitizeMol(mol)

    return mol


def readProtein(fname, **kwargs):
    mol = Chem.MolFromPDBBlock((gzip.open if fname.endswith('.gz') else open)(fname, 'rt').read(), **kwargs)
    # Chem.AssignAtomChiralTagsFromStructure(mol)
    # Chem.AssignStereochemistry(mol)
    # Chem.SanitizeMol(mol)

    return mol

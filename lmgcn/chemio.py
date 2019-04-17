import gzip
from rdkit import Chem

def readLigand(fname):
    return Chem.MolFromMolBlock((gzip.open if fname.endswith('.gz') else open)(fname, 'rt').read())

def readProtein(fname):
    return Chem.MolFromPDBBlock((gzip.open if fname.endswith('.gz') else open)(fname, 'rt').read())

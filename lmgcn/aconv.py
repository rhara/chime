"""
Atomic Conv for Complex inspired by Pande
"""

from .chemio import readLigand, readProtein
from .matrix import getAtomicNumVector, getDistanceMatrix
import sys
import numpy as np

__VERSION__ = '0.1.1'
__AUTHOR__ = 'hara@lifematics.co.jp'

np.set_printoptions(precision=3, threshold=np.inf, linewidth=300)

class ComplexFeaturizer:
    def __init__(self,
                 atomtype_list=[6,7,8,9,11,12,15,16,17,20,25,30,35,53],
                 radials=(1.6, 12.0, 0.5),
                 r_s=1.5,
                 sigma_s=1.0,
                 M=12,
                ):
        self.atomtype_list = atomtype_list
        self.radials = np.arange(radials[0], radials[1]+1e-7, radials[2])
        self.r_s = r_s
        self.sigma_s = sigma_s
        self.M = 12
        self.lig_fname = None
        self.pro_fname = None
        self.R = None
        self.Z = None
        self.E = None
        self.P = None

    def __call__(self, lig_fname, pro_fname):
        self.lig_fname = lig_fname
        self.pro_fname = pro_fname
        lig = readLigand(lig_fname)
        pro = readProtein(pro_fname)

        print('Ligand: %s: %d' % (lig_fname, lig.GetNumAtoms()), file=sys.stderr)
        print('Protein: %s: %d' % (pro_fname, pro.GetNumAtoms()), file=sys.stderr)

        self.neighborListConstruction(lig, pro)
        self.atomtypeConvolution()
        self.radialPooling()

    def neighborListConstruction(self, lig, pro):
        """
        Make R (distance) and Z (atomtype) of neighboring atoms
        Note: neighbor size is defined as self.M
        """

        N1 = lig.GetNumAtoms()
        N2 = pro.GetNumAtoms()

        atomtype1 = getAtomicNumVector(lig)
        atomtype2 = getAtomicNumVector(pro)

        self.R = np.zeros([N1, self.M])

        dist_mat = getDistanceMatrix(lig, pro)

        nbors = np.zeros([N1, self.M], dtype=np.int16)
        self.Z = np.zeros([N1, self.M], dtype=np.int16)
        for i in range(N1):
            nbors[i] = np.argsort(dist_mat[i])[:self.M]
            self.R[i] = dist_mat[i, nbors[i]]
            self.Z[i] = atomtype2[nbors[i]]

    def atomtypeConvolution(self):
        """
        Do the atom type convolution
        R(N,M), Z(N,M) -> E(N,M,Nat)
        """

        Na = len(self.atomtype_list)
        stacks = []
        for k in range(Na):
            atomtype = self.atomtype_list[k]
            Ka = (self.Z == atomtype).astype(np.int)
            v = self.R*Ka
            stacks.append(v)
        self.E = np.dstack(stacks)

    def radialPooling(self):
        """
        Do the radialPooling
        E(N,M,Nat) -> P(N,Nat,Nr)
        """

        def fc(r, cutoff):
            return (r < cutoff)*(np.cos(np.pi*r/cutoff) + 1)

        def fs(r, cutoff):
            return np.exp(-(r-r_s)*(r-r_s)*fc(r, cutoff)/sigma_s/sigma_s)

        beta_init = 1.0
        bias_init = 0.0

        r_s = self.r_s
        sigma_s = self.sigma_s
        beta = np.ones(len(self.radials))*beta_init
        bias = np.ones(len(self.radials))*bias_init

        N = self.E.shape[0]
        Na = len(self.atomtype_list)
        Nr = len(self.radials)

        self.P = np.zeros([N, Na, Nr])
        for i in range(Nr):
            cutoff = self.radials[i]
            r0_list = np.zeros(Na)
            v0 = beta[i]*np.sum(fs(r0_list, cutoff)) + bias[i]
            for j in range(N):
                for k in range(Na):
                    r_list = self.E[j,:,k]
                    if np.all(r_list == 0.):
                        self.P[j, k, i] = v0
                    else:
                        self.P[j, k, i] = beta[i]*np.sum(fs(r_list, cutoff)) + bias[i]

    def save(self, fname):
        """
        Save attributes in npz format
        """

        np.savez(fname,
                 lig_fname=self.lig_fname,
                 pro_fname=self.pro_fname,
                 R=self.R,
                 Z=self.Z,
                 E=self.E,
                 P=self.P)
        print('Saved in %s' % (fname,), file=sys.stderr)


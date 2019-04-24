"""
Atomic Conv for Complex inspired by Pande
"""

from .graph_features import getAtomTypeVector, getAtomPairDistanceMatrix, getNeighborList
import sys
import numpy as np

__VERSION__ = '0.1.1'
__AUTHOR__ = 'hara@lifematics.co.jp'


class ComplexFeaturizer:
    def __init__(self,
                 atomtype_list=[6,7,8,9,11,12,15,16,17,20,25,30,35,53],
                 radials=(1.6, 12.0, 0.5),
                 r_s=1.5,
                 sigma_s=1.0,
                 M=12,
                 padding=70,
                ):
        self.atomtype_list = atomtype_list
        self.radials = np.arange(radials[0], radials[1]+1e-7, radials[2])
        self.r_s = r_s
        self.sigma_s = sigma_s
        self.M = 12
        self.padding = padding
        self.lig = None
        self.pro = None
        self.R = None
        self.Z = None
        self.E = None
        self.P = None

    def __call__(self, lig, pro):
        self.lig = lig
        self.pro = pro
        lig_fname = lig.GetProp('fname')
        pro_fname = pro.GetProp('fname')
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

        atomtype1 = getAtomTypeVector(lig)
        atomtype2 = getAtomTypeVector(pro)

        self.R = np.zeros([N1, self.M])

        dist_mat = getAtomPairDistanceMatrix(lig, pro)
        self.R, self.Z = getNeighborList(lig, pro, M=self.M)

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

        if N < self.padding:
            _P = np.ones((self.padding - N, Na, Nr))*v0
            self.P = np.vstack((self.P, _P))

    def save(self, fname):
        """
        Save attributes in npz format
        """

        lig_fname = self.lig.GetProp('fname')
        pro_fname = self.pro.GetProp('fname')

        np.savez(fname,
                 lig_fname=lig_fname,
                 pro_fname=pro_fname,
                 R=self.R,
                 Z=self.Z,
                 E=self.E,
                 P=self.P)
        print('Saved in %s' % (fname,), file=sys.stderr)

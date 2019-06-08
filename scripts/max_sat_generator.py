import numpy as np
import matplotlib.pyplot as plt
import random as rm
import numpy.random as nrm

N = 5  # Number of bits
M = 20  # Number of clauses
K = 2   # Max-K-sat

# Generate a max-sat instance. Randomly select k-tuples of bits, assign
# a sign to each bit (to indicate negation), and create a clause from
# the k bits. Do this M times.
terms = []  # H = A_i*z_i + B_ij*z_i*z_j

for i in range(M):
    clause = rm.sample(range(N),K)
    signs  = 2*nrm.choice(range(2),K)-1
    terms.append((clause,signs))

def maxksat(z):
    H = len(terms)
    for i in range(len(terms)):
        clause=terms[i][0]
        signs =terms[i][1]
        p = 1
        for j in range(len(clause)):
            p = p*(1-signs[j]*z[clause[j]])
        H = H-p
    return H

print maxksat(nrm.choice(2,N))

import itertools
Nbit = list(itertools.product([0,1],repeat=N))
Elandscape = map(maxksat,Nbit)
plt.plot(Elandscape); plt.show()


# coding: utf-8

# In[72]:

import numpy as np
import numpy.random as nrm
import random as rm
from matplotlib import pyplot as plt
import time
import itertools


# File input

# In[405]:

clauses = []
f=open('/home/aniruddha/Documents/Projects/QuICS/AdOpt/ms_random/abrame-habet/max2sat/120v/s2v120c1200-1.cnf','r') # Input file, cnf format
while True:
    line = f.readline()
    if not line: break
    spline = line.split()
    if spline[0]=='c':
        continue
    elif spline[0]=='p':
        N = int(spline[2])
        M = int(spline[3])
    else:
        for c in spline[:-1]:
            if abs(int(c))>N or int(c)==0:
                print("Error: variable indices must be non-zero and less than %d"%N)
                break;
        if spline[-1]!='0':
            print("Error: clause descriptions must have a terminal 0")
            break;
        clauses.append(map(int,spline[:-1]))
if len(clauses)!=M: print("Error: Need %d clauses"%M)
f.close()


# Max sat function

# In[406]:

def maxksat(z):
    H = len(clauses)
    for i in range(len(clauses)):
        clause=map(abs,clauses[i])
        signs=map(np.sign,clauses[i])
        p = 1
        for j in range(len(clause)):
            p = p*(0.5*(1+signs[j])-signs[j]*z[clause[j]-1])
        H = H-p
    return H


# In[419]:

maxksat(nrm.choice(2,N))


# Exhaustive search to find maximum satisfiability (WARNING: SCALES EXPONENTIALLY IN # OF BITS)

# Nbit = list(itertools.product([0,1],repeat=N))
# minE = M
# minz = []
# for z in Nbit:
#     zE = maxksat(z)
#     if minE >= zE:
#         minE = zE
#         minz = z
# print minz, minE        

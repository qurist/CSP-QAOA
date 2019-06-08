#----------------------------------------------------------------------------#
# Title: A python helper library for QAOA-type optimization
# Author: Aniruddha Bapat
# Date: 05-28-2018
#
# Description: Here, I will maintain a library of frequently used functions
# in QAOA-type optimizations.
# Citation: Please cite this code if you use it in your research.
#----------------------------------------------------------------------------#

from ctypes import *
from FnLibrary import *
import scipy.sparse as sp
from scipy.optimize import basinhopping, minimize
import itertools as it
from numpy import pi

# Load c functions in qaoa
# I think this is OS-dependent. Create the appropriate shared c object
# for your OS (.so, .dll etc.) and load it into qc

scripts = '/home/aniruddha/Documents/Projects/QuICS/AdOpt/scripts'
qc      = cdll.LoadLibrary('%s/VQE/lib2local-qaoa.so'%scripts)
csp     = cdll.LoadLibrary('%s/VQE/libcsp.so'%scripts)

# States and Hamiltonian structs, ctyped into python
class state(Structure):
    _fields_ = [('n', c_int),
                ('N',c_int),
                ('realcur',POINTER(c_double)),
                ('imagcur',POINTER(c_double)),
                ('realbuf',POINTER(c_double)),
                ('imagbuf',POINTER(c_double))]
class ham(Structure):
    _fields_ = [('n',c_int),
                ('N',c_int),
                ('zzc',POINTER(c_double)),
                ('xc',POINTER(c_double)),
                ('zc',POINTER(c_double))]

# Bound class for basinhopping
class MaxTbound(object):
    """random displacement with bounds"""
    def __init__(self, T, stepsize=0.5):
        self.maxT= T
        

    def __call__(self, **kwargs):
        """take a random step but ensure the new position is within the bounds"""
        x = kwargs["x_new"]
        return (sum(x)<self.maxT)

    
# Set input and return types
qc.expectH.restype  = c_double
qc.overlap.resype   = c_double
qc.energy.restype   = c_double
qc.qaoa1energy.restype = c_double

csp.expectH.restype  = c_double
csp.energy.restype   = c_double

# Initialize state and hamiltonian
def initialize(n):
     psi     = state()
     H       = ham()
     qc.allocate(byref(psi),n)
     qc.allocateH(byref(H),n)
     return (psi, H)

# Generate Ham object from coefficients zzc, xc and zc (given as np arrays)
def hamGen(H, zzc, xc, zc):
     n = len(xc)
     H.zzc = (c_double * n**2)(*zzc.flatten().tolist())
     H.xc  = (c_double * n)(*xc.tolist())
     H.zc  = (c_double * n)(*zc.tolist())


# Function to compute Ham's k smallest eigenvalues and eigenvectors
def ground(H, Neigs):
     n = H.n
     energies= map(lambda i:qc.energy(i,byref(H)), range(1<<n))
     Ham     = np.diag(energies)
     for i in range(1<<n):
          for j in range(n):
               Ham[i,i^(1<<j)] = H.xc[j]
     return sp.linalg.eigsh(Ham, k=Neigs, which='SA')

# Inner product magnitude of two states given as complex np arrays
def overlap(psi1, psi2): 
     return np.abs(np.vdot(psi1,psi2))

# Convert C state object into numpy array of computational basis amplitudes
def c2pyState(psi):
     n = psi.n
     return np.array(psi.realcur[:1<<n]) + 1.0j*np.array(psi.imagcur[:1<<n])
 
# Overlap of psi with the Neigs smallest ground states of H
def gsOverlap(psi, H, Neigs):
     val, vec = ground(H, Neigs)
     final    = c2pyState(psi)
     retvalue = 0
     for k in range(Neigs):
          retvalue += overlap(final, vec[:,k])**2
     return np.sqrt(retvalue)

#betagamma = all betas first, then all gammas
def expectQAOAp(psi, H, betagamma):
    ppsi = pointer(psi)
    pH   = pointer(H)
    csp.uniform(ppsi)
    p = len(betagamma)/2
    csp.evolveQAOA(ppsi, pH, (c_double * p)(*betagamma[:p]), (c_double * p)(*betagamma[p:]), p)
    return csp.expectH(psi, pH)

# Perform an inductive local optimization of QAOA angles. The code returns
# optimal angles and energy. Psi is now the result of the optimal evolution.
def optQAOA(psi, H, pmax, typeOfOpt='BFGS'):
    ppsi = pointer(psi)
    pH   = pointer(H)
    fOpt = lambda bg: expectQAOAp(psi, H, bg.tolist())
    bg0  = 0.5*np.ones(2)
    bgCur= bg0   
    for p in range(1,pmax+1):
        bgNew  = np.concatenate((bgCur[:p-1], [bgCur[p-1]], bgCur[p-1:], [bgCur[-1]]))
        opt    = minimize(fOpt, bg0 if p==1 else bgNew, method=typeOfOpt)
        bgCur  = opt.x
        E      = expectQAOAp(psi, H, bgCur.tolist())
        #if E!=opt.fun: return "Error: Energy expectation does not match optimized value"
    return (bgCur, E) 

def optQAOAglobal(psi, H, p, Nit=10, Temp=1.0, Tmax=20):
    ppsi = pointer(psi)
    pH   = pointer(H)
    fOpt = lambda bg: expectQAOAp(psi, H, bg.tolist())
    bg0  = 2.0*np.ones(2*p)
    bounds = [(0,pi/2)]*p + [(0,np.inf)]*p
    minKwargs = dict(method="L-BFGS-B", bounds=bounds)
    stepCond  = MaxTbound(Tmax)
    opt  = basinhopping(fOpt, bg0, niter=Nit, T=Temp, minimizer_kwargs=minKwargs, accept_test=stepCond)
    bgOpt= opt.x
    E      = expectQAOAp(psi, H, bgOpt.tolist())
    return (opt.x, E)

# Full search of angles, with resolution of M for every pi/2 angles
# beta range: [0,pi/2]
# Gamma range: [0, 3*pi]

def optQAOAgreedy(psi, H, pmax, typeOfOpt='BFGS'):
    ppsi = pointer(psi)
    pH   = pointer(H)
    bg0  = [0.5]*2
    bgOpt  = []
    for p in range(pmax):
        fOpt   = lambda bg: expectQAOAp(psi, H, np.concatenate((bgOpt,bg)).tolist())
        opt    = minimize(fOpt, bg0, method=typeOfOpt)
        bgOpt  += opt.x.tolist()
    E    = expectQAOAp(psi, H, bgOpt)
    return (np.array(bgOpt), E) 

def optQAOAfull(psi, H, p, M=10, trunc=-np.inf):
    S    = 2 # scale factor gamma:beta
    bMax = pi/2
    gMax = S*bMax
    indOpt= np.inf
    Opt  = np.inf
    
    configs = it.product(*([np.linspace(0,bMax, M,endpoint=False)]*p + [np.linspace(0,gMax, S*M,endpoint=False)]*p))

    for ind, i in enumerate(configs):
        Ecur = expectQAOAp(psi, H, i)
        if trunc!=np.inf and Ecur < trunc:
            Opt = Ecur
            indOpt = ind
            break
        indOpt= indOpt if Opt < Ecur else ind       
        Opt  = min(Opt, Ecur)
    configs = it.product(*([np.linspace(0,bMax, M,endpoint=False)]*p + [np.linspace(0,gMax, 6*M,endpoint=False)]*p))
    return (list(configs)[indOpt], Opt)

def printOut(psi, val, vec, bgOpt, Eopt):
    final = c2pyState(psi)
    final2gs = 0
    sumval   = 0
    Neigs    = len(vec.T)
    for k in range(Neigs):
        sumval += overlap(final, vec[:,k])**2
    final2gs=np.sqrt(sumval)

    final2initial = 1./psi.N*overlap(final, np.ones(psi.N))
    for i in range(len(bgOpt)):
        print "%5.4f "%bgOpt[i],
    print
    #print "%5.4f %5.4f"%(Eopt, val[0])
    print "%7.6f %5.4f"%(np.abs(Eopt/(val[0])), final2gs)
    print "%5.4f"%final2initial


# Symmetrize a state w.r.t to the global Z2 symmetry present in 2-local Hamiltonians
# Only physical when you are in the symmetric sector
def z2symmetrize(psi):
     sym = psi + psi[::-1]
     return sym/np.sqrt(overlap(sym,sym))

# Entanglement entropy
# Input: a state psi (as an np.arrar) in the computational basis, and the cut site k
def entent(psi, k):
     N = len(psi)
     n = len(format(N,'b'))-1 # effectively n = log2(N) for powers of two
     m = n-k
     prho = np.zeros((1<<k,1<<k),dtype=np.complex_) # partial density matrix
     x = np.arange(1<<k)
     for i in range(1<<m):
          ppsi = psi[x*(1<<m)+i] # Partial state vector
          prho += np.outer(ppsi,np.conj(ppsi))
     val, vec = eigh(prho)
     return sum(map(lambda x: -x*np.log(abs(x)+0.000000000001), val))

# The analytical formula for energy under QAOA1
# betagamma = beta first, then gamma
def expectQAOA1(H, betagamma): 
    return qc.qaoa1energy(byref(H), c_double(betagamma[0]), c_double(betagamma[1]))



##############################################################################
# ZZ + X + Z Hamiltonian coefficient generators
# Output format: a list of lists of the form ((Z, ZZ, ..), (X, XX, ..))
# So, for a 2 local Hamiltonian, expect ((Z, ZZ), (X,))

# The standard 2-local long-range model with power law interactions
def lr(n, alpha, J, B):
    Ki  = np.zeros(n)
    Jij = np.array([[(J/(abs(i-j))**alpha if i!=j else 0) for i in range(n)] for j in range(n)])
    Bi  = B*np.ones(n)
    return ((Ki, Jij), (Bi,))

# Read parameters and prepare coefficients
def readParams():
    return

# Construct an unweighted MaxSat instance from a cnf file
# NOTE1: cnf Clauses must be 1-indexed.
# NOTE2: This method is pretty memory-inefficient for
# large clause sizes. Ideal for up to ~4Sat

def readCnf(path2file):
    f=open(path2file,'r')
    clauses   = []
    maxLength = 0
    while True:
        line = f.readline()
        if not line: break
        spline = line.split()
        if spline[0]=='c':
            continue
        elif spline[0]=='p':
            n = int(spline[2])
            M = int(spline[3])
        else:
            for c in spline[:-1]:
                if abs(int(c))>n or int(c)==0:
                    print("Error: variable indices must lie between 0 and %d"%n)
                    break;
            if spline[-1]!='0':
                print("Error: clause descriptions must have a terminal 0")
                break;
            maxLength = max(maxLength, len(spline)-1)
            clauses.append(map(int,spline[:-1]))
    if len(clauses)!=M: print("Error: Need %d clauses"%M)
    f.close()
    const = 0
    Bi    = np.zeros(n)
    Jall  = [np.zeros([n]*i) for i in range(1,maxLength+1)]
    for i in range(M):
        clause=map(lambda x: abs(x)-1,clauses[i])
        signs=map(np.sign,clauses[i])
        K = len(clause)
        for k in range(1, K+1):
            for inds in it.combinations(range(K), k):
                Jall[k-1][inds] += (1./(1<<K))*np.prod([signs[x] for x in inds])
        const += 1/(1<<K)
             
    return (Jall, (Bi,))
        

    

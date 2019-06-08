from ctypes import *
from FnLibrary import *
from qaoaLibrary import *
import scipy.sparse as sp
from scipy.optimize import basinhopping, minimize
import itertools as it
from numpy import pi
import sys
import time
from matplotlib.pyplot import *

# Load c functions in qaoa
# I think this is OS-dependent. Create the appropriate shared c object
# for your OS (.so, .dll etc.) and load it into qc

start = time.time()

scripts = '/home/aniruddha/Documents/Projects/QuICS/AdOpt/scripts'
csp      = cdll.LoadLibrary('%s/VQE/libcsp.so'%scripts)

# States and Hamiltonian structs, ctyped into python
class state(Structure):
    _fields_ = [('n', c_int),
                ('N',c_int),
                ('realcur',POINTER(c_double)),
                ('imagcur',POINTER(c_double)),
                ('realbuf',POINTER(c_double)),
                ('imagbuf',POINTER(c_double))]
class cspham(Structure):
    _fields_ = [('n',c_int),
                ('N',c_int),
                ('m',c_int),
                ('maxsat',c_double),
                ('clauses',POINTER(POINTER(c_int))),
                ('lengths',POINTER(c_int)),
                ('weights',POINTER(c_double)),
                ('diag', POINTER(c_double))]
H = cspham()

csp.expectH.restype  = c_double
csp.energy.restype   = c_double

csp.readcnf(sys.argv[1],byref(H))

psi  = state()
ppsi = pointer(psi)
pH   = pointer(H)

csp.allocate(ppsi,H.n)
csp.uniform(ppsi)
#csp.energize(pH)
#print "Average energy: ", csp.expectH(psi, pH)


Emin, sat = ground_csp(H)


energies = map(lambda x: csp.energy(x, byref(H)), np.arange(psi.N))

# Isolate states with a certain cost
costc = np.arange(psi.N)[np.array(energies) == 30]

x = 5
c = energies[x]
for i in range(5):
    nbhd = []
    x = costc[i]
    for inds in it.combinations(range(psi.n), 2):
        nbhd.append(energies[(x^(1<<inds[0]))^(1<<inds[1])])
#        nbhd.append(energies[x^(1<<inds[0])])

    hist(nbhd);show()
    
    

print time.time()-start

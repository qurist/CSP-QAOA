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
#print H.m, H.n, H.maxsat
#for i in range(H.m):
#    for j in range(H.lengths[i]):
#        print H.clauses[i][j],
#    print 
#print csp.energy(int(sys.argv[2]), byref(H))    

psi  = state()
ppsi = pointer(psi)
pH   = pointer(H)

csp.allocate(ppsi,H.n)
csp.uniform(ppsi)
csp.energize(pH)
print "Average energy: ", csp.expectH(psi, pH)

p=int(sys.argv[2])
bg   = [pi/2-pi/15]*p + [0.4*pi, 0.4*pi]*(p/2)
csp.evolveQAOA(ppsi, pH,(c_double * p)(*bg[:p]),(c_double * p)(*bg[p:]),p)


Emin, sat = ground_csp(H)
bg, Eopt  = optQAOA_csp(psi, H, p, sat)

#plot(bg[:p]); plot(bg[p:]);show()
print bg
final = abs(c2pyState(psi))**2
perm  = final

print final[sat]
print sum(final[sat])

for i in range(1<<psi.n):
    perm[i] = final[int(((psi.n-len(bin(i)[2:]))*'0' + bin(i)[2:])[::-1],2)]

#print ground_csp(H)

#plot(final)
energies = map(lambda x: csp.energy(x, byref(H)), range(psi.N))
#plot(energies)
#show()

pbycost = np.zeros(H.m)

for i in range(psi.N):
    pbycost[int(energies[i])]+= final[i]

    
#plot(pbycost)
#xlim(0,20)
#hist(final, bins=energies[] s=0.1)
#show()

print time.time()-start

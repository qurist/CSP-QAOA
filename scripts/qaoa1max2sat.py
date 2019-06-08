from ctypes import *
from FnLibrary import *
from qaoaLibrary import *
import sys
from scipy.optimize import minimize, basinhopping, differential_evolution
import time
from numpy import pi

start = time.time()


# Read input description of 2SAT formula
nvars       = sys.argv[1]#'120'
nclauses    = sys.argv[2]
index       = sys.argv[3]

#path = '/home/aniruddha/Documents/Projects/QuICS/AdOpt/ms_random/abrame-habet/max2sat/' + nvars + 'v/s2v' + nvars + 'c' + nclauses + '-' + index + '.cnf'

path = '/home/aniruddha/Documents/Projects/QuICS/AdOpt/MaxSAT-instances/ms_crafted/maxcut/dimacs-mod/hamming10-2.clq.cnf'

#path = '/home/aniruddha/Documents/Projects/QuICS/AdOpt/MaxSAT-instances/instance/urand_' + nvars + '_' + nclauses + '.wcnf'

#path = '/home/aniruddha/Documents/Projects/QuICS/AdOpt/MaxSAT-instances/instance/bal-'+ nvars + '-' + nclauses + '.wcnf'


#############################################################################
# MAXSAT FORMULA INPUT

f=open(path,'r')

clauses = []
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
                print("Error: variable indices cannot exceed %d"%N)
                break;
        if spline[-1]!='0':
            print("Error: clause descriptions must have a terminal 0")
            break;
        clauses.append(map(int,spline[:-1]))
if len(clauses)!=M: print("Error: Need %d clauses"%M)
f.close()

#############################################################################
# CLAUSES ARE 1-INDEXED, WHILE THE QUBITS ARE 0-INDEXED!!!

psi, H  = initialize(n)

zzc = np.zeros(n**2)
zc  = np.zeros(n)
xc  = np.zeros(n)

for i in range(M):
    clause=map(abs,clauses[i])
    signs=map(np.sign,clauses[i])
    J = 1
    # Remember, clauses are 1-indexed!
    for j in range(2):
        zc[clause[j]-1] += signs[j]
        J*=signs[j]
    zzc[(clause[0]-1)*n+clause[1]-1]+= J
    zzc[(clause[1]-1)*n+clause[0]-1]+= J
s = [1.0]*n

hamGen(H, zzc, zc, xc)

fOpt = lambda bg: expectQAOA1(H, bg)

opt     = minimize(fOpt, [0.3,0.2], method='BFGS')
#opt     = basinhopping(fOpt, [0.3,0.2], niter=5, T=3.0)
#opt      = differential_evolution(fOpt, bounds)

print nclauses + '-' + index, opt.x[0]%(pi/2), opt.x[1]%(pi), 0.25*(M+opt.fun), 0.75 - opt.fun/(4*M), time.time() - start

#qc.deallocate(ppsi)
#qc.deallocateH(ph)

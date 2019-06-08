import pandas as pd
from pandas.plotting import scatter_matrix
import matplotlib.pyplot as plt
from sklearn import model_selection
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC

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


csp.expectH.restype  = c_double
csp.energy.restype   = c_double

H = cspham()
psi  = state()
ppsi = pointer(psi)
pH   = pointer(H)

NINSTANCES = 1000#int(sys.argv[1])
p = 1#int(sys.argv[2])

path = 'Random3Sat/'

### INSTANCE DICTIONARY. CONTAINS ESSENTIAL DESCRIPTORS
instances = dict({'Nvars'        : [],
                  'Nclauses'     : [],
                  'Cdensity'     : [],
                  'IncidenceMean': [],
                  'IncidenceStd' : [],
                  'Unsat'        : []})

for i in range(1,p+1):
    instances['gamma'+'%d'%i]= []
    instances['beta'+'%d'%i] = []

### READ INSTANCES
for inst in range(NINSTANCES):
    csp.readcnf(path + '%d.cnf'%inst, byref(H))
    csp.allocate(ppsi,H.n)
    csp.energize(pH)
    csp.uniform(ppsi)
    
    instances['Nvars'].append(H.n)
    instances['Nclauses'].append(H.m)
    instances['Cdensity'].append(H.m/float(H.n))
    PosIncidence= np.zeros(H.n) # number of times a literal appears
    NegIncidence= np.zeros(H.n) # number of times a negated literal appears
    PosIncidence= np.zeros(H.n) # number of times a literal appears
    Incidence= np.zeros(H.n) # total number of times a literal appears
    
    for i in range(H.m):
        for j in range(H.lengths[i]):
            a = H.clauses[i][j]
            if a>0:
                PosIncidence[a-1] += 1
            elif a<0:
                NegIncidence[-a-1] += 1
    Incidence = PosIncidence + NegIncidence
    PosIncidenceMean = np.mean(PosIncidence)
    NegIncidenceMean = np.mean(NegIncidence)
    PosIncidenceStd = np.std(PosIncidence)
    NegIncidenceStd = np.std(NegIncidence)
    IncidenceMean = np.mean(Incidence)
    IncidenceStd = np.std(Incidence)

    instances['IncidenceMean'].append(IncidenceMean)
    instances['IncidenceStd'].append(IncidenceStd)
    
    Emin, sat = ground_csp(H)
    BGopt, Eopt = optQAOA_csp(psi, H, p, sat)
    for i in range(p):
        instances['beta' +'%d'%(i+1)].append(BGopt[i])
        instances['gamma'+'%d'%(i+1)].append(BGopt[p+i])
    instances['Unsat'].append((Eopt - Emin)/float(H.m))
    csp.deallocate(ppsi)
    csp.deallocateH(pH)


### NOW, WE LEARN
qaoaset = pd.DataFrame(instances)
print(qaoaset.shape)
    



print time.time()-start


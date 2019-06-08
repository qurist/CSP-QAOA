from qaoaLibrary import *
import sys
import time
from matplotlib import pyplot as plt
from numpy import random as nrm
from scipy.interpolate import interp1d


#### MAIN FUNCTION, WILL BE CALLED BELOW ##############
def main(filepath, p, OptMethod):
    # Initialize state and hamiltonian
    psi = state()
    H   = cspham()
    csp.readcnf(filepath, byref(H))
    csp.allocate(byref(psi), H.n)
    csp.energize(byref(H))
    ##### QAOA PATH  #################################
    opt    = minimize(lambda bg: expectQAOAp_csp(psi, H, bg.tolist()), 0.1*np.ones(2), method=OptMethod)
    bg0    = np.abs(opt.x)
    E0     = expectQAOAp_csp(psi, H, opt.x.tolist())
    if p==1: return (bg0, E0)
    
    #### p>1, use path interpolation #################
    # Optimize energy
    fOpt = lambda bg: expectQAOAp_csp(psi, H, bg.tolist())

    # Initialize path function
    Bpath = lambda i: bg0[0]-i*0.2
    Gpath = lambda i: bg0[1]+i*0.2
    q0=20

    for q in range(2,p+1):
        spacing = np.linspace(0,1,q)
        bgCur = np.array([Bpath(i) for i in spacing] + \
                         [Gpath(i) for i in spacing])
        opt   = minimize(fOpt, bgCur, method=OptMethod)
        bgCur = np.copy(opt.x)
        
        Bpath = interp1d(spacing,bgCur[:q], kind='cubic' if q>3 else 'linear')
        Gpath = interp1d(spacing,bgCur[q:], kind='cubic' if q>3 else 'linear')
        if q>q0:
            plt.plot(spacing, bgCur[:q],'r', alpha=((q-q0)/(float(p)-q0))); plt.plot(spacing, bgCur[q:],'b', alpha=((q-q0)/(float(p)-q0)))


    # Return
    plt.show()
    #plt.plot(spacing, bgCur[:p]); plt.plot(spacing, bgCur[p:]); plt.show()
    Eopt = expectQAOAp_csp(psi, H, bgCur.tolist())
    return (opt.x, Eopt)

########## BEGIN CODE ################
start_time = time.time()

# Read in parameters
filepath  = sys.argv[1]
p  = int(sys.argv[2])

OptMethod = 'Powell'
bg, E = main(filepath,p,OptMethod)

#### INSTEAD, YOU COULD JUST GIVE AN INITIAL GUESS ###
#### FOR THE PATH AND FEED IT INTO AN OPTIMIZER  #####
# bg = np.array([0.48039303,  0.32333747,  0.26701939,  0.23536333,  0.21655129,  0.20081751,
#              0.19619766,  0.1903008,   0.18237125,  0.18187484,  0.17746348,  0.17415601,
#              0.16898964,  0.1606949,   0.08454086,  0.08484914,  0.1933862,   0.23924315,
#              0.28206927,  0.30033203,  0.31720761,  0.32813725,  0.34102611,  0.35333936,
#              0.36759798,  0.37426245,  0.38132785,  0.39121731,  0.39758966,  0.37712917])
# opt = minimize(lambda bg: expectQAOAp(psi, H, bg.tolist()), bg, method=OptMethod)
# Eopt = expectQAOAp(psi, H, opt.x.tolist())

#### WRITE OUTPUT ####
print bg
print E
# bgfile = open('QAOApath_p=%d.dat'%p,'a')
# for i in range(len(bg)):
#     bgfile.write("%6.5f "%bg[i])
# bgfile.write('\n')
# bgfile.close()
# Efile  = open('QAOAenergy_p=%d.dat'%p, 'a')
# Efile.write("%f\n"%E)
# Efile.close()
print(time.time() - start_time)

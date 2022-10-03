import numpy as np

def L1norm(C, Cref, dr, R):
    return 1/R*np.sum(dr*abs((C-Cref)))

def L2norm(C, Cref, dr, R):
    return np.sqrt(1/R*np.sum((dr*(C-Cref)**2)))

def Linfnorm(C, Cref):
    return np.max(abs(C-Cref))

import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.path.join("src"))

# Modules du programme
from solver import *
from analyticalSolution import *

if __name__ == "__main__":

    # Inputs
    n = 10
    Deff = 10e-10
    S = 10e-8
    k = 4e-9
    R = 0.5
    Ce = 10
    dt = 1.0

    dr = R/(n-1.0)
    # Stockage des paramètres afin de faciliter l'accès
    parameters = [n, dr, dt, Deff, S, k, Ce]

    # Initial condition
    C0 = np.zeros(n)

    # Résolution
    # C = unsteadySolve(parameters, 1000, C0, constantSource=True, secondOrder=True)
    C = steadySolve(parameters, constantSource=True, secondOrder=True)
    print(C)

    # Solution analytique
    Ca = analyticalSolution(S, Deff, R, Ce, n)
    print(Ca)

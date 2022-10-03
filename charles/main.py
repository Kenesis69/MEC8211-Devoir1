import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.path.join("src"))

# Modules du programme
from solver import *
from analyticalSolution import *
from error import *

if __name__ == "__main__":

    # base inputs
    n = 5
    Deff = 10e-10
    S = 10e-8
    k = 4e-9
    R = 0.5
    Ce = 10
    dt = 1.0

    dr = R/(n-1.0)
    # Stockage des paramètres afin de faciliter l'accès
    parameters = [n, dr, dt, Deff, S, k, Ce]

    # Condition initiale
    C0 = np.zeros(n)

    ############################
    # Analyse de convergence en
    # régime stationnaire et
    # terme source constant
    ############################

    # Analyse de convergence - ordre 2
    with open("convergence_ordre2.dat", "w") as file:
        file.write("CONVERGENCE DU SHÉMA D'ORDRE 2 EN ESPACE\n")
        file.write("NPOINTS L1 L2 LINF\n")
        for i in range(0,10):
            # Calcul des solutions
            C = steadySolve(parameters)
            Cref = analyticalSolution(S, Deff, R, Ce, n)

            # Calcul des erreurs avec la solution analytique
            L1 = L1norm(C, Cref, dr, R)
            L2 = L2norm(C, Cref, dr, R)
            Linf = Linfnorm(C, Cref)

            file.write("%i %.12f %.12f %.12f\n" % (n, L1, L2, Linf))

            n*=2
            dr = R/(n-1.0)
            parameters[0] = n
            parameters[1] = dr

    # Le nombre de point est remis à 5 pour la prochaine analyse
    n = 5
    parameters[0] = n

    # Analyse de convergence - ordre 1
    with open("convergence_ordre1.dat", "w") as file:
        file.write("CONVERGENCE DU SHÉMA D'ORDRE 1 EN ESPACE\n")
        file.write("NPOINTS L1 L2 LINF\n")
        for i in range(0,10):
            # Calcul des solutions
            C = steadySolve(parameters, secondOrder = False)
            Cref = analyticalSolution(S, Deff, R, Ce, n)

            # Calcul des erreurs avec la solution analytique
            L1 = L1norm(C, Cref, dr, R)
            L2 = L2norm(C, Cref, dr, R)
            Linf = Linfnorm(C, Cref)

            file.write("%i %.12f %.12f %.12f\n" % (n, L1, L2, Linf))

            n*=2
            dr = R/(n-1.0)
            parameters[0] = n
            parameters[1] = dr

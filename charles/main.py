import numpy as np
import matplotlib.pyplot as plt

# Calcul de la matrice A du système linéaire Ax=b
#       *Note... il y a trois "switch" (valeur 0 ou 1) en entrée qui permettent
#               d'utiliser la même fonction pour résoudre plusieurs cas possibles
#               (i.e stationnaire ou instationnaire, premier ou second ordre,
#               terme source constant ou du premier ordre)
def calcA(parameters, steadySwitch=0, constantSourceSwitch=0, secondOrderSwitch=0):
    n = parameters[0]
    dr = parameters[1]
    dt = parameters[2]
    Deff = parameters[3]
    S = parameters[4]
    k = parameters[5]

    A = np.zeros((n,n))

    # Condition limite à r=0 (Neumann)
    A[0,0] = -1.0
    A[0,1] = 1.0

    # Centre du domaine
    r = 0
    for i in range(1,n-1):
        r += dr
        A[i,i-1] = -r*dt*Deff + secondOrderSwitch*0.5*dr*dt*Deff
        A[i,i] = steadySwitch*r*dr**2 + (1-secondOrderSwitch)*dr*dt*Deff + 2*r*dt*Deff + \
            constantSourceSwitch*k*dt*dr**2*r
        A[i,i+1] = -(-0.5*secondOrderSwitch+1)*dr*dt*Deff - r*dt*Deff

    # Condition limite à r=R (Dirichlet)
    A[n-1,n-1] = 1.0

    # print(A)
    return A

# Calcul du vecteur b du système linéaire Ax=b
#       *Note... il y a deux "switch" (valeur 0 ou 1) en entrée qui permettent
#               d'utiliser la même fonction pour résoudre plusieurs cas possibles
#               (i.e stationnaire ou instationnaire, terme source constant ou du premier ordre)
def calcRHS(parameters, Cold, steadySwitch=0, constantSourceSwitch=0):
    n = parameters[0]
    dr = parameters[1]
    dt = parameters[2]
    Deff = parameters[3]
    S = parameters[4]
    k = parameters[5]

    b = np.zeros(n)

    # Condition limite à r=0 (Neumann)
    b[0] = 0.0

    # Centre du domaine
    r = 0
    for i in range(1,n-1):
        r += dr
        b[i] = steadySwitch*Cold[i]*r*dr**2 - (1-constantSourceSwitch)*S*dr**2*dt*r

    # Condition limite à r=R (Dirichlet)
    b[n-1] = Ce

    # print(b)
    return b

# Solveur instationnaire contenant la boucle temporelle
def unsteadySolve(parameters, tend, C0, constantSource=True, secondOrder=True):

    # "booléen" qui permait de d'activer les termes instationnaires
    steadySwitch = 1 # Si régime instationnaire

    # "booléen" qui permait d'activer/désactiver des termes selon le terme source
    if (constantSource):
        constantSourceSwitch = 0 # Si source constante
    else:
        constantSourceSwitch = 1 # Si source est une réaction de premier ordre

    # "booléen" qui permait d'activer/désactiver les termes de deuxième ordre
    if (secondOrder):
        secondOrderSwitch = 1
    else:
        secondOrderSwitch = 0

    dt = parameters[2]
    t = 0
    Cold = C0
    while t<tend:
        A = calcA(parameters, steadySwitch, constantSourceSwitch, secondOrderSwitch)
        b = calcRHS(parameters, Cold, steadySwitch, constantSourceSwitch)
        Cold = np.linalg.solve(A, b)
        t+=dt
        print(Cold)

    return Cold

# Solveur stationnaire résolvant le système linéaire en un seul coup
def steadySolve(parameters, constantSource=True, secondOrder=True):

    # "booléen" qui permait de désactiver les termes instationnaires
    steadySwitch = 0 # Si régime stationnaire
    parameters[2] = 1.0 # On impose un dt de 1 pour annuler l'effet sur les termes stationnaires

    # "booléen" qui permait d'activer/désactiver des termes selon le terme source
    if (constantSource):
        constantSourceSwitch = 0 # Si source constante
    else:
        constantSourceSwitch = 1 # Si source est une réaction de premier ordre

    # "booléen" qui permait d'activer/désactiver les termes de deuxième ordre
    if (secondOrder):
        secondOrderSwitch = 1
    else:
        secondOrderSwitch = 0

    n = parameters[0]
    A = calcA(parameters, steadySwitch, constantSourceSwitch, secondOrderSwitch)
    b = calcRHS(parameters, np.zeros(n), steadySwitch, constantSourceSwitch)
    C = np.linalg.solve(A, b)

    return C

# Solution analytique sur un domaine donné
def analyticalSolution(S, Deff, R, Ce, n):
    r = np.linspace(0, R, n)
    return (0.25*S/Deff*R**2*(r**2/R**2-1) + Ce)


#############################################################
#############################################################


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
    parameters = [n, dr, dt, Deff, S, k]

    # Initial condition
    C0 = np.zeros(n)

    # Résolution
    # C = unsteadySolve(parameters, 1000, C0, constantSource=True, secondOrder=True)
    C = steadySolve(parameters, constantSource=True, secondOrder=True)
    print(C)

    # Solution analytique
    Ca = analyticalSolution(S, Deff, R, Ce, n)
    print(Ca)

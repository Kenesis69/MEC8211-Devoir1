import numpy as np

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
    # S = parameters[4]
    k = parameters[5]
    # Ce = parameters[6]

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
    # Deff = parameters[3]
    S = parameters[4]
    # k = parameters[5]
    Ce = parameters[6]

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

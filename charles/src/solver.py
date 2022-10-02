import numpy as np
from linearSystem import *

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

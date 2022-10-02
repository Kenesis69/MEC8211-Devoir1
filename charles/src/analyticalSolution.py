import numpy as np

# Solution analytique sur un domaine donné
def analyticalSolution(S, Deff, R, Ce, n):
    r = np.linspace(0, R, n)
    return (0.25*S/Deff*R**2*(r**2/R**2-1) + Ce)

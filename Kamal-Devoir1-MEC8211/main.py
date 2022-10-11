from matplotlib.pyplot import plot
from Fonction import *

import matplotlib.pyplot as plt

class prm:
    D_eff = 10**-10
    k=4*10**-9
    S_anal = 0.29
    
    C_i = 10


dt = 100000000
t = 1000000000000
n = 5
R = 0.5
a , b = diff_stationnaire_degre1(prm,100,R)


a , c = diff_stationnaire_degre2(prm,100,R)


rep = diff_transitoire_degre1(prm,n,R,dt,t)
print(rep)

rep = diff_transitoire_degre2(prm,n,R,dt,t)
print(rep)
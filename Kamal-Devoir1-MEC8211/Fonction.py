#fonction avec le terme source a S = 10**-8

import numpy as np

    #D_eff = 10**-10
    #k=4*10**-9
    #S_anal = 10**-8
# CECI EST UNE MODIFICATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

def diff_stationnaire_degre1(prm,n,R):
    r = np.linspace(0,0.5,n)
    dr = R/(n-1)
    A = np.zeros([n,n])
    b = np.linspace(0,0,n)
    for i in range(0,n):
        alpha = prm.D_eff/dr**2
        beta = -prm.D_eff/(r[i]*dr) - 2*prm.D_eff/dr**2 - prm.k
        ceta = prm.D_eff/(r[i]*dr) + prm.D_eff/dr**2
        if i == n-1:
            A[i,i] = 1
            b[i] = prm.C_i
        elif i == 0:
            A[i,i] = -1
            A[i,i+1] = 1
            b[i] = 0
        else: 
            A[i,i-1] = alpha
            A[i,i] = beta
            A[i,i+1] = ceta
            b[i] = 0
    Rep = np.linalg.solve(A,b)
    return r , Rep


def diff_stationnaire_degre2(prm,n,R):
    r = np.linspace(0,0.5,n)
    dr = R/(n-1)
    A = np.zeros([n,n])
    b = np.linspace(0,0,n)
    for i in range(0,n):
        alpha = prm.D_eff/dr**2 -prm.D_eff/(r[i]*dr) 
        beta = - 2*prm.D_eff/dr**2 - prm.k
        ceta = prm.D_eff/(r[i]*dr) + prm.D_eff/dr**2
        if i == n-1:
            A[i,i] = 1
            b[i] = prm.C_i
        elif i == 0:
            A[i,i] = -1
            A[i,i+1] = 1
            b[i] = 0
        else: 
            A[i,i-1] = alpha
            A[i,i] = beta
            A[i,i+1] = ceta
            b[i] = 0
    Rep = np.linalg.solve(A,b)
    return r , Rep


def diff_transitoire_degre1(prm,n,R,dt,t):
    A = np.zeros([n,n])
    b  = np.linspace(0,0,n)
    r = np.linspace(0,0.5,n)
    dr = R/(n-1)
    time = 0
    while time < t:
        for i in range(0,n):
            alpha = dt*prm.D_eff/dr**2
            beta = -dt*prm.D_eff/(r[i]*dr) -2*dt*prm.D_eff/dr**2 -dt*prm.k -1
            ceta = dt*prm.D_eff/(r[i]*dr) + dr*prm.D_eff/dr**2
            if i == n-1:
                A[i,i] = 1
                b[i] = prm.C_i
            elif i == 0:
                A[i,i] = -1
                A[i,i+1] = 1
                b[i] = 0
            else: 
                A[i,i-1] = alpha
                A[i,i] = beta
                A[i,i+1] = ceta
                b[i] = -b[i]
        Rep = np.linalg.solve(A,b)
        time = time + dt
    return Rep

def diff_transitoire_degre2(prm,n,R,dt,t):
    A = np.zeros([n,n])
    b  = np.linspace(0,0,n)
    r = np.linspace(0,0.5,n)
    dr = R/(n-1)
    time = 0
    while time < t:
        for i in range(0,n):
            alpha = dt*prm.D_eff/dr**2-dt*prm.D_eff/(r[i]*dr)
            beta =  -2*dt*prm.D_eff/dr**2 -dt*prm.k -1
            ceta = dt*prm.D_eff/(r[i]*dr) + dr*prm.D_eff/dr**2
            if i == n-1:
                A[i,i] = 1
                b[i] = prm.C_i
            elif i == 0:
                A[i,i] = -1
                A[i,i+1] = 1
                b[i] = 0
            else: 
                A[i,i-1] = alpha
                A[i,i] = beta
                A[i,i+1] = ceta
                b[i] = -b[i]
        Rep = np.linalg.solve(A,b)
        time = time + dt
    return Rep

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 13:08:25 2018

@author: bartlomiejkos
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 11:31:57 2018

@author: bartlomiejkos
"""
import numpy as np
import matplotlib.pyplot as plt
import cmath

me = 9.109383e-31 # kg  # declaration of variables
hbar = 1.053571e-34 # Js
eV = 1.602176e-19 # J
mw = me * 0.067
mb = me * 0.1
#E = 0.1
V = 0.334 * eV
b = 10 * 1e-9
#w = 10 * 1e-9
z = cmath.sqrt(-1)
Ener = []
Tran = []
w = 9 * 1e-9
for e in range(31400, 31450):   # iterating over specified range of energies
    if e == 33400:
        continue
    else:
        e = e * 0.00001
        k = cmath.sqrt(2*mw*e*eV/hbar**2)
        K = cmath.sqrt(2*mb*(V-e*eV)/hbar**2)
        M1 = np.ones((2,2), dtype='complex')
        M1[1][0] = cmath.sqrt(-1) * k / mw
        M1[1][1] = -cmath.sqrt(-1) * k / mw
        M2 = np.ones((2,2), dtype='complex')
        M2[1][0] = K / mb
        M2[1][1] = -K / mb
        M3 = np.ones((2,2), dtype='complex')
        M3[0][0] = cmath.exp(K * b)
        M3[0][1] = cmath.exp(-K * b)
        M3[1][0] = K * cmath.exp(K * b) / mb
        M3[1][1] = -K * cmath.exp(-K * b) / mb
        M4 = np.ones((2,2), dtype='complex')
        M4[0][0] = cmath.exp(cmath.sqrt(-1) * k * b)
        M4[0][1] = cmath.exp(-cmath.sqrt(-1) * k * b)
        M4[1][0] = cmath.sqrt(-1) * k * cmath.exp(cmath.sqrt(-1) * k * b) / mw
        M4[1][1] = -cmath.sqrt(-1) * k * cmath.exp(-cmath.sqrt(-1) * k * b) / mw
        M5 = np.ones((2,2), dtype='complex')
        M5[0][0] = cmath.exp(cmath.sqrt(-1) * k * (b + w))
        M5[0][1] = cmath.exp(-cmath.sqrt(-1) * k * (b + w))
        M5[1][0] = cmath.sqrt(-1) * k * cmath.exp(cmath.sqrt(-1) * k * (b + w)) / mw
        M5[1][1] = -cmath.sqrt(-1) * k * cmath.exp(-cmath.sqrt(-1) * k * (b + w)) / mw
        M6 = np.ones((2,2), dtype='complex')
        M6[0][0] = cmath.exp(K * (b + w))
        M6[0][1] = cmath.exp(-K * (b + w))
        M6[1][0] = K * cmath.exp(K * (b + w)) / mb
        M6[1][1] = -K * cmath.exp(-K * (b + w)) / mb
        M7 = np.ones((2,2), dtype='complex')
        M7[0][0] = cmath.exp(K * (2 * b + w)) 
        M7[0][1] = cmath.exp(-K * (2 * b + w))
        M7[1][0] = K * cmath.exp(K * (2 * b + w)) / mb
        M7[1][1] = -K * cmath.exp(-K * (2 * b + w)) / mb
        M8 = np.ones((2,2), dtype='complex')
        M8[0][0] = cmath.exp(cmath.sqrt(-1) * k * (2 * b + w))
        M8[0][1] = 0
        M8[1][0] = cmath.sqrt(-1) * k * cmath.exp(cmath.sqrt(-1) * k * (2 * b + w)) / mw
        M8[1][1] = 0
        M = np.zeros((2,2), dtype='complex')
        M = np.matmul(np.linalg.inv(M1), M2)
        M = np.matmul(M, np.linalg.inv(M3))
        M = np.matmul(M, M4)
        M = np.matmul(M, np.linalg.inv(M5))
        M = np.matmul(M, M6)
        M = np.matmul(M, np.linalg.inv(M7))
        M = np.matmul(M, M8)
        T = 1 / (M[0][0] * np.conj(M[0][0]))
        T = float(T)
        Tran.append(T)
        Ener.append(e)
plt.plot(Ener, Tran, label="w=%s nm"%(round(w*1e9,2),))
Tran = [] 
Ener = []
plt.legend(loc='lower right', fontsize='x-small')
plt.yscale('log')
plt.show()


# sklejanie funkcji falowej

E = 0.27745
m0 = 9.10938356
mw = 0.067 * m0
mb = 0.10 * m0
V0 = 0.334
hbar = 1.0545718

k = cmath.sqrt(2*mw*E/hbar**2)
K = cmath.sqrt((2*mb*(V0-E))/hbar**2)
def wave_free(F, R, k, x):
    return F * np.exp(cmath.sqrt(-1) * k * x) + R * np.exp(-cmath.sqrt(-1) * k * x)
def wave_in(F, R, K, x):
    return F * np.exp(K * x) + R * np.exp(-K * x)



MKL = np.ones((2,1), dtype='complex')
K1 = 1
L = 0
MKL[0] = K1
MKL[1] = L
MAB = np.ones((2,1), dtype='complex')
MAB = np.matmul(M, MKL)
A = MAB[0]
B = MAB[1]
MCD = np.ones((2,1), dtype='complex')
MFG = np.ones((2,1), dtype='complex')
MHJ = np.ones((2,1), dtype='complex')
MHJ = np.matmul(np.linalg.inv(M7), M8)
MHJ = np.matmul(MHJ, MKL)
H = MHJ[0]
J = MHJ[1]
MFG = np.matmul(np.linalg.inv(M5), M6)
MFG = np.matmul(MFG, MHJ)
F = MFG[0]
G = MFG[1]
MCD = np.matmul(np.linalg.inv(M3), M4)
MCD = np.matmul(MCD, MFG)
C = MCD[0]
D = MCD[1]
print(A,B,C,D,F,G,H,J,K1,L)
a = 5
b = 10
w = 10
z1=8
z2=12
z3=25
psi = []
z = np.linspace(-15, 40, 1000)

for zi in z:
    if zi < 0:
        F = wave_free(A, B, k, zi)
        f = np.real(F * np.conj(F))
        psi.append(f)
    elif 0 < zi < z1:
        F = wave_in(C, D, K, zi)
        f = np.real(F * np.conj(F))
        psi.append(f)
    elif z1 < zi < z2:
        F = wave_free(F, G, k, zi)
        f = np.real(F * np.conj(F))
        psi.append(f)
    elif z2 < zi < z3:
        F = wave_in(H, J, K, zi)
        f = np.real(F * np.conj(F))
        psi.append(f)
    elif z3 < zi:
        F = wave_free(K1, L, k, zi)
        f = np.real(F * np.conj(F))
        psi.append(f)
S = []
for item in psi:
    s = float(item)
    S.append(s)

plt.plot(z, S)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 11:31:57 2018

@author: bartlomiejkos
"""
import numpy as np
import pandas as pd
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
for w in range(2, 11):
    w = w * 1e-9
    for e in range(25000, 36000):   # iterating over specified range of energies
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
    plt.plot(Ener, Tran, label="w=%s"%(w*1e9,))
    Tran = [] 
    Ener = []
plt.legend(loc='lower right', fontsize='x-small')
plt.yscale('log')
plt.show()

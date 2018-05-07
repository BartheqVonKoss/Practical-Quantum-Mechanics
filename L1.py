#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 24 16:54:15 2018

@author: bartlomiejkos
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



'''---------Applied Quantum Mechanics---------'''

# Zadanie1
'''
Zbadać zależność czterech pierwszych poziomów energetycznych od szerokości studni w dla 
elektronu w nieskończonej studni potencjału w GaAs (m*=0.067m). Rozważyć szerokości 
1nm <w <30nm. Wyniki przedstawić na wykresie.
'''

def zadanie1():
    me = 9.109383e-31 # kg  # declaration of variables
    hbar = 1.053571e-34 # Js
    eV = 1.602176e-19 # J
    n = 4
    mx = me * 0.067
    E = np.empty([n, 29])  # definition of 2D matrix to hold results
    for ni in range(1, n + 1):   # iterating over energy lvls
        for wi in range(-14, 15):  # iterating over well's width
            if wi == 0:
                continue
            else:
                Ei = (hbar**2 * np.pi**2 * ni**2) / (8 * mx * (wi/2)**2 * eV * 1e-9**2)  # i-th energy lvl
                E[ni-1][wi+15-1] = Ei
    E = np.transpose(E)
    E=pd.DataFrame(E, columns=('n1','n2','n3','n4'))  # wrangling data
    plt.plot(range(-14,15), E['n1'], 'ro')
    plt.plot(range(-14,15), E['n2'], 'g+')
    plt.plot(range(-14,15), E['n3'], 'y^')
    plt.plot(range(-14,15), E['n4'], 'bs')
    plt.xlabel('width (nm)')
    plt.ylabel('energy (eV)')
    plt.xlim(-18,18)
    plt.show()
# Zadanie2

'''
Skończona studnia kwantowa szerokości o 10nm . Znaleźć numerycznie wszystkie poziomy 
energetyczne dla głębokości V0=0.5 eV, 1eV i 2eV. Przyjąć masę ładunku m*=0.067m.    
'''

def zadanie2():
    me = 9.109383e-31 # kg  # declaration of variables
    hbar = 1.053571e-34 # Js
    eV = 1.602176e-19 # J
    V = [0.5, 1 , 2] # potential lvls
    E = np.zeros([8,3])
    j = 0
    for i in V: # iterating over potential lvls
        a = 1
        ni = 1
        while a == 1:
            Ei = (hbar**2 * np.pi**2 * ni**2) / (2 * me * 0.067 * 10**2 * eV * 1e-9**2)  # i-th energy lvl
            if Ei < i:    # checking whether Ei is still smaller than V0
                #print(Ei, i)    # energy lvl in accord to V0
                E[ni-1][j] = Ei
                ni = ni + 1
            else:
                a = 0
        j = j + 1
    PLOT(E)
def PLOT(E):
    plt.plot([0, 2], [E[0][0], E[0][0]], color='r', linestyle='--', linewidth=2)
    plt.plot([0, 2], [E[1][0], E[1][0]], color='r', linestyle='--', linewidth=2)
    plt.plot([2, 4], [E[0][1], E[0][1]], color='c', linestyle=':', linewidth=2)
    plt.plot([2, 4], [E[1][1], E[1][1]], color='c', linestyle=':', linewidth=2)
    plt.plot([2, 4], [E[2][1], E[2][1]], color='c', linestyle=':', linewidth=2)
    plt.plot([2, 4], [E[3][1], E[3][1]], color='c', linestyle=':', linewidth=2)
    plt.plot([4, 6], [E[0][2], E[0][2]], color='b', linestyle='-.', linewidth=2)
    plt.plot([4, 6], [E[1][2], E[1][2]], color='b', linestyle='-.', linewidth=2)
    plt.plot([4, 6], [E[2][2], E[2][2]], color='b', linestyle='-.', linewidth=2)
    plt.plot([4, 6], [E[3][2], E[3][2]], color='b', linestyle='-.', linewidth=2)
    plt.plot([4, 6], [E[4][2], E[4][2]], color='b', linestyle='-.', linewidth=2)
    plt.plot([0, 6], [0, 0], color='k', linestyle='-', linewidth=3)
    plt.ylabel('energy (eV)')
    plt.show()

# Zadanie3
'''
Obliczyć współczynnik transmisji T przez barierę energetyczną dla następujących danych 
mw=0.067 m0, mb=0.10 m0, V0= 0.334 eV. Zbadać T(E) dla kilku szerokości bariery w,
 2nm< w<10nm. Wyniki obliczeń przedstawić na wykresach.
'''

def zadanie3():
    #    eV = 1.602176e-19 # J  # declaration of variables
    m0 = 9.10938356
    mw = 0.067 * m0
    mb = 0.10 * m0
    V0 = 0.334
    hbar = 1.0545718
    E = []
    w = [3, 5, 8]
    
    for wi in w:
        T = []
        E = []
        for e in np.linspace(0.01, 1.6, 10000):
            if e == V0:
                continue
            elif e < V0:
                E.append(e)
                k1 = np.sqrt((2 * mw * e) / (hbar**2))
                k2 = np.sqrt((2 * mb * (V0 - e)) / hbar**2)
                t = 1 / (1 + 0.25 * (k1 / k2 + k2 / k1)**2 * (np.sinh(k2 * wi))**2)
                T.append(t)
            else:
                E.append(e)
                k1 = np.sqrt((2 * mw * e) / (hbar**2))
                k2 = np.sqrt((2 * mb * (e - V0)) / hbar**2)
                t = 1 / (1 + 0.25 * (k1 / k2 - k2 / k1)**2 * (np.sin(k2 * wi))**2)
                T.append(t)
        plt.plot(E, T, label = "w =%s nm" %wi)
    plt.legend()
    plt.grid()
    plt.xlabel('E [eV]')
    plt.ylabel('Transmission coefficient \nT [0-1]')
    plt.show()


#zadanie1()
#zadanie2()    
#zadanie3()

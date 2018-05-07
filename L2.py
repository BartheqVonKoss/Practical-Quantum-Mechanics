import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.optimize import newton

'''-------Applied Quantum Mechanics-------'''

# Zadanie1

'''
Zastosować zmodyfikowaną metodę strzałów do znalezienia energii i funkcji falowych
'''

def Zadanie1():
    n = 1000
    w = 3
    z = np.linspace(-w, w, n)
    dz = 2 * w / n
    
    E = np.linspace(0, 50, 1000)
    V = np.zeros(n)
    V = 33.4*V
    
    m = np.ones(n)
    mw = m * 0.067
    
    psi = np.zeros(n)
    chi = np.ones(n)

    ans = []
    eps = 1e-1
    _E = 0
    for Ei in E:
        wf = PSI(dz, z, psi, chi, V, mw, Ei)
        if wf[-1] * wf[-2] < 0:
            dE = Ei - _E
            if dE > eps:
                ans.append((Ei, wf.copy()))
    
    
    for item in ans:
        plt.plot(z, item[1], label='E=%.3f' % item[0])
    
    plt.grid()
    plt.title("Wave functions for infinite potential well")
    plt.legend()
    plt.show()
def PSI(dz, z, psi, chi, V, m, E):
    for i in range(1, len(z)):
        psi[i] = psi[i-1] + dz*m[i]*chi[i-1]
        chi[i] = chi[i-1] + dz * 2.0 * (V[i] - E) * psi[i-1]
    return psi



# Zadanie2

'''
Zastosować zmodyfikowaną metodę strzałów do znalezienia energii i funkcji falowych:
(a) dla pojedynczej studni o szerokości w
(b) dwóch studni rozdzielonych barierą o szerokości b i potencjale V0. 
Przyjąć mw=0.067 m0, mb=0.10 m0, V0= 0.334 eV.
'''
def Zadanie2():
    m0 = 9.10938356
    mw = 0.067 * m0
    mb = 0.10 * m0
    V0 = 0.334
    hbar = 1.0545718
    ws = [3,4,5,6]
    b = 3
    N = 200
    
    def f(E, alpha):
        k = np.sqrt(2 * mw * E / hbar**2)
        K = np.sqrt(2 * mb * (V0 - E) / hbar**2)
        A = np.cos(w * k) * np.cosh(b * K)
        fact = (mb**2 * k**2 - mw**2 * K**2) / (2 * mb * mw * k * K)
        B = np.sin(w * k) * np.sinh(b * K)
        return A - fact * B - np.cos(alpha * 2 * np.pi)
    
    x = np.linspace(-1, 1, N)
    for w in ws: 
        e = []
        for alpha in x:
            s = newton(f, 0.08, args=(alpha,))
            e.append(s)
        plt.plot(x, e, label="w=%s nm"%w)
    plt.legend(loc='upper center', fontsize='small')
    plt.title('Dependency of energy vs wave vector g')
    plt.xlabel('g')
    plt.ylabel('E (eV)')
    plt.grid()
    plt.show()
    
    alpha = 1
    ene = np.linspace(0.001, 0.3, N)
    z = []
    for E in ene:
        c = f(E, alpha)
        z.append(c)
    plt.plot(ene, z)
    plt.title('f(E) with g = 1')
    plt.xlabel('E (eV)')
    plt.axhline(linewidth=1, color='black')
    plt.grid()
    plt.show()
    
    
#Zadanie1()
#Zadanie2()

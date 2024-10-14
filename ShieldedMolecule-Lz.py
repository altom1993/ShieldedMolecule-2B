# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 09:57:00 2024

@author: yang
"""

import numpy as np
from scipy.integrate import simps
from scipy.optimize import root_scalar
import matplotlib 
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


# constant parameters
eps0 = 8.854187817*10**(-12) # Vacuum permittivity：farad/m
debye = 3.3356409519815204*10**(-30) # Debye, coulomb*m
hbar = 1.054571817*10**(-34) # Plank constant hbar, the unit is Js
Hz = 1 # 1/s

# give the parameters of molecules
d = 475/100 # dipolar moment of NaCs, with the unit of Debye
mu =2.58870105526947*10**(-25)/2 # reduced mass of two NaCs molecules, the unit is kg
# d = 285/100 # dipolar moment of NaK, with the unit of Debye
# mu = 1.0619412082166079*10**(-25)/2 # reduced mass of two NaK molecules, the unit is kg

# give the value of detuning(delta_test) and trap frequency(omega_test), the unit is Hz
delta_test = 10*2*np.pi*10**6*Hz
omega_test = 1*2*np.pi*10**6*Hz
# give the range of the ratio between Rabi frequency and detuning, Or donates Omega/Delta
rangeOr = [1/2, 10]



#%%
# dimensionless
length = (debye**2/(4*np.pi*eps0*hbar*Hz))**(1/3)

def l0(delta):
    return (d**2/delta)**(1/3)*length
def E0(delta):
    return hbar**2/(2*mu*(l0(delta))**2)
def Delta0(delta):
    return hbar*delta/E0(delta)
def aomega(omega):
    return np.sqrt(hbar/(mu*omega))
def aomega0(delta, omega):
    return aomega(omega)/l0(delta)

# express the internal couplings 
def h2rho3(omegar,rho):
    a = omegar**2
    b = 1 + a
    c = np.sqrt(b)
    d = 1 + 1/c
    e = 1 - 1/c
    f = rho**3
    return np.array([[-a/(12*b), 0, d*np.sqrt(e)/4, -omegar/(6*np.sqrt(2)*b), 0,-np.sqrt(d)*e/4, a/(12*b)],[0, 1/6*(d-3*(1+c)*f), 0, 0, -omegar/(6*c), 0, 0], [d*np.sqrt(e)/4, 0, 1/12*(-d-6*(1+c)*f), np.sqrt(d/2)/(2*c), 0, omegar/(12*c), -d*np.sqrt(e)/4], [-omegar/(6*np.sqrt(2)*b), 0, np.sqrt(d/2)/(2*c), -1/(6*b)-c*f, 0, -np.sqrt(e/2)/(2*c), omegar/(6*np.sqrt(2)*b)], [0, -omegar/(6*c), 0, 0, 1/6*(e-3*(f+3*c*f)), 0, 0], [-np.sqrt(d)*e/4, 0, omegar/(12*c), -np.sqrt(e/2)/(2*c), 0, 1/12*(-e-6*(f+3*c*f)), np.sqrt(d)*e/4], [a/(12*b), 0, -d*np.sqrt(e)/4, omegar/(6*np.sqrt(2)*b), 0, np.sqrt(d)*e/4, 1/12*(-1+1/b-24*c*f)]])



## define the microwave-shielded potential
def vbo(omegar, delta, rho):
    value, vector= np.linalg.eigh(h2rho3(omegar,rho))
    return Delta0(delta)*1/rho**3*np.max(value)
## define the gauge potential
operator_a = np.diag([0, -1, -2, 0, -1, -2, 0])
def aphi(omegar, rho):
    value, vector = np.linalg.eigh(h2rho3(omegar,rho))
    index = np.argmax(value)
    vector1 = vector[:,index]
    return 1/rho*(vector1@operator_a@vector1)
## express the harmonic trap potential
def vtrap(delta,omega,rho):
    return rho**2/aomega0(delta, omega)**4
## express the total potenial
def veff(m, omegar, delta, omega, rho):
    return (m/rho-aphi(omegar, rho))**2 - 1/(4*rho**2) + vtrap(delta, omega, rho) + vbo(omegar, delta, rho)



#%%
from pack.potential1D import Potential1D

def normallizeV(V: callable, w, N, min=0.5, energy_cut=0):
    """
    cut the divergent part of the potential

    right cut: if we need to cut the right part of the potential.
    """

    T = 2*np.pi/w
    print('*******************', min)
    if energy_cut == 0:
        energy_cut = N*w*10
    print(energy_cut, '========================')
    r0 = root_scalar(lambda x: V(x)-energy_cut, method='bisect', x0=1e-10,
                     bracket=[1e-1, min]).root
    r1 = np.inf

    def normedV(x):
        x += 1e-10
        get_outer = -np.sign((x-r0)*(r1-x))/2 + 1/2
        get_inner = 1 - get_outer
        return get_outer*energy_cut + get_inner*V(x)
    points = []
    for r in [r0, r1]:
        if r < T:
            points.append(r)
    return normedV, points




# calculate the relative angular momentum for the ground state for different Omega/Delta
listOr = np.linspace(rangeOr[0], rangeOr[1], 20)
listE0 = np.zeros(len(listOr))
listLz = np.zeros(len(listOr))

for i in range(0,len(listOr)):
    Or = listOr[i]
    
    w = 2 # determine the cutoff of rho
    N = 601 # N determine the cutoff of the basis
    levels = 1 # number of bound state to be calculated

    def vtotal0(rho): # total potenial for m=0
        return veff(0, Or, delta_test, omega_test, rho)
    def vtotaln1(rho): # total potenial for m=-1
        return veff(-1, Or, delta_test, omega_test, rho)

    nV0, points0 = normallizeV(lambda r: vtotal0(rho=r), w=w, N=N)
    nVn1, pointsn1 = normallizeV(lambda r: vtotaln1(rho=r), w=w, N=N)

    p3d0 = Potential1D(V=lambda x:nV0(x), N=N, w=w, verbose=1, points=points0, zero='left')# Numercically solve the bound state and the eigen wave functions for relative motion in the potential vtotal0.
    p3dn1 = Potential1D(V=lambda x:nVn1(x), N=N, w=w, verbose=1, points=pointsn1, zero='left')# Numercically solve the bound state and the eigen wave functions for relative motion in the potential vtotaln1.
    
    vals0 = p3d0.get_eigenvals(levels=levels) # energy of the lowest bound state for m=0
    valsn1 = p3dn1.get_eigenvals(levels=levels) # energy of the lowest bound state for m=-1
    
    listrho = np.linspace(0, 2*np.pi/w, 1000)
    listaphi_0 = np.array([listrho[i]*aphi(Or, listrho[i]) for i in range(1, len(listrho))]) # list of A_phi
    listaphi = np.insert(listaphi_0, 0, 0) # A_phi(rho) is singular at rho=0. here, we simply take rho*A(rho)=0 at rho=0, and this affect barely the result of Lz, since the wave funtion at rho=0 is zero.

    if vals0 < valsn1: # determine which is the bound state, m=0 or -1.
        wf = p3d0.get_eigenwf()
        Lz = -1*simps(wf[0]*wf[0]*listaphi, listrho) # calculate the relative angular momentum for the ground state.
    else:
        wf = p3dn1.get_eigenwf()
        Lz = -1 - 1*simps(wf[0]*wf[0]*listaphi, listrho)
    
    listLz[i] = Lz

    print(f"Omega_r = {Or}")
    print(f"normalization = {simps(wf[0]*wf[0], listrho)}") #check normalization
    print(f"Lz = {Lz}")

print(listLz)


#%%
np.savetxt('./Lz.txt',listLz, fmt='%.8f')











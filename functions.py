from scipy.optimize import fsolve
import numpy as np
from classes_hx import water


def f(Re, eD):
    A = eD/3.7+(6.7/Re)**0.9
    fo = 1/(-2*np.log10(eD/3.7-5.02/Re*np.log10(A)))**2

    if eD:
        f = fsolve(lambda x: 1/x**0.5+2.0*np.log10(eD/3.7+2.51/Re/x**0.5), fo)
    else:
        f = fsolve(lambda x: 1/x**0.5-2.0*np.log10(Re*x**0.5)+0.8, fo)
    return f[0]

def V(m_dot,liquid,area):
    rho = liquid.rho
    V = m_dot/(rho*A)
    return V

def Re(V,d,mu,rho):
    Re = (V*d*rho)/mu
    return Re

def hi(V_tube,di,liquid):

    mu = liquid.mu
    rho = liquid.rho
    Pr = liquid.Pr
    k = liquid.k

    Re = Re(V_tube,di,mu,rho)

    Nu = 0.023*Re**0.8*Pr**0.3

    hi = Nu*k/di

    return hi

def ho(V_shell,do,liquid,tube_layout='t'):
    if tube_layout == 't':
        c = 0.2
    elif tube_layout == 's':
        c = 0.15
    else:
        c= 0.2
        print('error, invalid tube layout')

    mu = liquid.mu
    rho = liquid.rho
    Pr = liquid.Pr
    k = liquid.k

    Re = Re(V_shell,do,mu,rho)

    Nu = c*Re**0.6*Pr**0.3

    ho = Nu*k/do

    return ho


def LMTD(T_1in,T_2in,T_1out,T_2out):
    lmtd = ((T_2in - T_1out) - (T_2out - T_1in)) / np.log((T_2in - T_1out) / (T_2out - T_1in))
    return lmtd

def U_inside(hi,ho,Ai,Ao,ri,ro,k_copper,L):

    U = 1 / (1/hi + (Ai*np.log(ro/ri))/(2*np.pi*k_copper*L) + Ai/(Ao*ho))
    return U


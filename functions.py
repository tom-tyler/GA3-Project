from math import log10
 
from scipy.optimize import fsolve

def f_chen(Re, eD):
    A = eD/3.7+(6.7/Re)**0.9
    f = 1/(-2*log10(eD/3.7-5.02/Re*log10(A)))**2
    return f

def f_colebrook(Re, eD):
    fo = f_chen(Re, eD)
    if eD:
        f = fsolve(lambda x: 1/x**0.5+2.0*log10(eD/3.7+2.51/Re/x**0.5), fo)
    else:
        f = fsolve(lambda x: 1/x**0.5-2.0*log10(Re*x**0.5)+0.8, fo)
    return f[0]
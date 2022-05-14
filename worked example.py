import numpy as np
from iapws._iapws import _Liquid
from classes_hx import HX

from functions import f_colebrook

baffles = 9
length = 0.35
tubes = 13

k_tube = 386
baffle_spacing = length/(baffles+1)

hot_water=_Liquid(273.15+40, 0.101325)            #Ice at normal melting point
# print(hot_water["rho"])
# print(hot_water["k"])
# print(hot_water["mu"])
# print(hot_water["cp"])
# print(hot_water["cp"]*hot_water["mu"]/hot_water["k"])

hx = HX(13,9,14e-3,350e-3,500e-3,16e-4)

print(hx.A_shell)
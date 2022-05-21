from hx_functions import mdot_dP
from hx_classes import water

h_w = water(60,46)

m = mdot_dP(0.06e5,'h',h_w,2022)
#print(m)
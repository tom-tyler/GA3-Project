import numpy as np
from iapws._iapws import _Liquid
from hx_classes import HX,water

import hx_functions as hxf

baffles = 9
length = 0.35
tubes = 13

k_tube = 386



hx = HX(13,9,14e-3,350e-3,500e-3,16e-4,tube_layout='s',shell_passes=1)

#initial values:
m_h = 0.45 #initial guess for hot mass flow rate
m_c = 0.50 #initial guess for cold mass flow rate

#specified inlet temperatures
T_inh = 60
T_inc = 20
#initial guesses for outlet temperatures
T_outh = 50
T_outc = 25

#initial guesses for inlet and outlet pressures
P_inh = 1.01325e5
P_outh = 1.01325e5
P_inc = 1.01325e5
P_outc = 1.01325e5

#creating hot and cold water objects
h_w = water(T_inh,T_outh,P_inh,P_outh)
c_w = water(T_inc,T_outc,P_inc,P_outc)

#mass flow in one tube
m_tube = m_h/hx.tube_number

#various velocities needed
V_tube = hxf.V(m_tube,h_w,hx.tube.c_area)
V_nozzle_h = hxf.V(m_h,h_w,hx.nozzle_c_area)
V_shell = hxf.V(m_c,c_w,hx.A_shell)
V_nozzle_c = hxf.V(m_c,c_w,hx.nozzle_c_area)

#sigma, used to find ke and kc
sigma = hx.sigma
# need to write code to give ke and kc still, so will add something here to do this

#pressure drop in a single tube
dP_tube = hxf.dP_tube(hx.tube_length,hx.tube.d_inner,h_w,V_tube)

#for now, letting kc = 0.45, and ke = 0.8 (correct for 13 tubes)
#pressure drop due to separation leaving tubes
dP_in_plus_out = hxf.dP_inout(h_w,V_tube,Kc = 0.45,Ke = 0.8)

#head loss in nozzle
dP_nozzles_h = 2 * hxf.dP_nozzle(V_nozzle_h,h_w)

#overall pressure drop
dP_tube_ovr = hx.tube_number*dP_tube + dP_in_plus_out + dP_nozzles_h


# now need iteration routine to get m_h such that dP_tube_ovr matches figure 6 from handout

#cold side
N = hx.shell_passes*hx.tube_number

dP_shell = hxf.dP_shell(V_shell,c_w,hx.tube.d_outer,N,tube_layout='s')

dP_nozzles_c = 2 * hxf.dP_nozzle(V_nozzle_c,c_w)

dP_shell_ovr = dP_shell + dP_nozzles_c

#now need iteration routine to get m_c such that dP_shell_ovr matches figure 6 from handout


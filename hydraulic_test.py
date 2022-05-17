import numpy as np
from iapws._iapws import _Liquid
from hx_classes import HX,water

import hx_functions as hxf


accuracy = 0.001
K_hmix = 6.0
K_baffles = 7.5

hx = HX(tube_number = 13,
        baffle_number = 14,
        pitch = 12e-3,
        tube_length = 362e-3,
        shell_length = 450e-3,
        baffle_area = 16e-4,
        tube_layout='t',
        shell_passes=1,
        nozzle_bore=25e-3)

#initial values:
m_h = 0.45 #initial guess for hot mass flow rate
m_c = 0.50 #initial guess for cold mass flow rate

#specified inlet temperatures
T_inh = 53.4
T_inc = 19.7

#initial guesses for outlet temperatures
T_outh = 46.8
T_outc = 26.1

#initial guesses for inlet and outlet pressures
P_inh = 1.01325e5
P_outh = 1.01325e5
P_inc = 1.01325e5
P_outc = 1.01325e5

#initial heat transfer parameters
heat_transfer,eff = 1,1
heat_transfer_new,eff_new = 1,1
Q_counter = 0
per_e_Q,per_e_eff = 1,1

#creating hot and cold water objects
h_w = water(T_inh,T_outh,P_inh,P_outh)
c_w = water(T_inc,T_outc,P_inc,P_outc)

m_counter = 0
dP_e_h = 1
m_e_h = 1
dP_e_c = 1
m_e_c = 1


while (abs(dP_e_h) > accuracy) and (abs(m_e_h) > accuracy) and (abs(dP_e_c) > accuracy) and (abs(m_e_c) > accuracy):

    #heat capacities
    Cc = m_c*c_w.cp
    Ch = m_h*h_w.cp

    #mass flow in one tube
    m_tube = m_h/hx.tube_number

    #various velocities needed
    V_tube = hxf.V(m_tube,h_w,hx.tube.c_area)
    V_nozzle_h = hxf.V(m_h,h_w,hx.nozzle_c_area)
    V_shell = hxf.V(m_c,c_w,hx.A_shell)
    V_nozzle_c = hxf.V(m_c,c_w,hx.nozzle_c_area)

    #sigma, used to find ke and kc
    sigma = hx.sigma
    sigma_nozzle = hx.sigma_nozzle
    print('sigma nozzle:',sigma_nozzle)

    #pressure drop in a single tube
    dP_tube = hxf.dP_tube(hx.tube_length,hx.tube.d_inner,h_w,V_tube)

    #pressure drop due to separation leaving tubes
    dP_in_plus_out = hxf.dP_inout(h_w,V_tube,sigma)

    #pressure drop due to separation leaving nozzle
    dP_in_plus_out_nozzle = hxf.dP_inout(h_w,V_nozzle_h,sigma_nozzle)

    #head loss in nozzle
    dP_nozzles_h = 2 * hxf.dP_nozzle(V_nozzle_h,h_w)

    dP_mixing = K_hmix*h_w.rho*((V_nozzle_h+V_tube)/2)**2

    #overall pressure drop
    dP_tube_ovr = dP_tube + dP_in_plus_out + dP_nozzles_h + dP_in_plus_out_nozzle + dP_mixing

    # now need iteration routine to get m_h such that dP_tube_ovr matches figure 6 from handout
    m_h, dP_e_h, m_e_h = hxf.mdot_dP(m_h,dP_tube_ovr,'h')
    
    print(f'm_h: {m_h},dP: {dP_tube_ovr}')

    #cold side
    N = hx.shell_passes*hx.tube_number

    dP_shell = K_baffles*hxf.dP_shell(V_shell,c_w,hx.tube.d_outer,N,tube_layout='s')

    dP_nozzles_c = 2 * hxf.dP_nozzle(V_nozzle_c,c_w)

    dP_shell_ovr = dP_shell + dP_nozzles_c

    #now need iteration routine to get m_c such that dP_shell_ovr matches figure 6 from handout
    m_c, dP_e_c, m_e_c = hxf.mdot_dP(m_c,dP_shell_ovr,'c')
    
    print(f'm_c: {m_c},dP: {dP_shell_ovr}')

    m_counter += 1
    if m_counter > 100:
            print('exceeded max iterations for m,dP')
            break

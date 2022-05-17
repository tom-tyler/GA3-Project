import numpy as np
from iapws._iapws import _Liquid
from hx_classes import HX,water

import hx_functions as hxf


accuracy = 0.001

hx = HX(13,14,12e-3,362e-3,450e-3,16e-4,tube_layout='t',shell_passes=1,nozzle_bore=25e-3)

#initial values:
m_h = 0.45 #initial guess for hot mass flow rate
m_c = 0.50 #initial guess for cold mass flow rate

#specified inlet temperatures
T_inh = 53.4
T_inc = 19.7

#initial guesses for outlet temperatures
T_outh = 50
T_outc = 25

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

while (abs(per_e_Q) > accuracy) and (abs(per_e_eff) > accuracy):

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
        print(sigma_nozzle)
    
        #pressure drop in a single tube
        dP_tube = hxf.dP_tube(hx.tube_length,hx.tube.d_inner,h_w,V_tube)

        #pressure drop due to separation leaving tubes
        dP_in_plus_out = hxf.dP_inout(h_w,V_tube,sigma)

        #pressure drop due to separation leaving nozzle
        dP_in_plus_out_nozzle = hxf.dP_inout(h_w,V_nozzle_h,sigma_nozzle)

        #head loss in nozzle
        dP_nozzles_h = 2 * hxf.dP_nozzle(V_nozzle_h,h_w)

        #overall pressure drop
        dP_tube_ovr = dP_tube + dP_in_plus_out + dP_nozzles_h + dP_in_plus_out_nozzle

        # now need iteration routine to get m_h such that dP_tube_ovr matches figure 6 from handout
        m_h, dP_e_h, m_e_h = hxf.mdot_dP(m_h,dP_tube_ovr,'h')
      
        print(f'm_h: {m_h},dP: {dP_tube_ovr}')

        #cold side
        N = hx.shell_passes*hx.tube_number

        dP_shell = hxf.dP_shell(V_shell,c_w,hx.tube.d_outer,N,tube_layout='s')

        dP_nozzles_c = 2 * hxf.dP_nozzle(V_nozzle_c,c_w)

        dP_shell_ovr = dP_shell + dP_nozzles_c

        #now need iteration routine to get m_c such that dP_shell_ovr matches figure 6 from handout
        m_c, dP_e_c, m_e_c = hxf.mdot_dP(m_c,dP_shell_ovr,'c')
      
        print(f'm_c: {m_c},dP: {dP_shell_ovr}')

        m_counter += 1
        if m_counter > 100:
                print('exceeded max iterations for m,dP')
                break


    #thermal design
    hi = hxf.hi(V_tube,hx.tube.d_inner,h_w)
    ho = hxf.ho(V_shell,hx.tube.d_outer,c_w,tube_layout='s')

    U = hxf.U_inside(hi,ho,hx.tube.d_inner,hx.tube.d_outer,hx.tube_length)

    A_con = hx.convection_area

    #now solve thermal design equations by iteration to get T_outh and T_outc. also find P_outh and P_outc

    T_outc_new,T_outh_new = 1,1
    T_counter = 0
    per_eh,per_ec = 1,1
    while (abs(per_ec) > accuracy) and (abs(per_eh) > accuracy):
        T_outc_new = T_inc + (1/Cc)*A_con*U*hxf.LMTD(T_inc,T_inh,T_outc,T_outh)
        T_outh_new = T_inh - (1/Ch)*A_con*U*hxf.LMTD(T_inc,T_inh,T_outc,T_outh)
        per_ec = (T_outc_new - T_outc)/T_outc
        per_eh = (T_outh_new - T_outh)/T_outh
        T_outc = T_outc_new
        T_outh = T_outh_new
        T_counter += 1

        if T_counter > 100:
            print('exceeded max iterations for T')
            break

    heat_transfer_new = hxf.Q(Cc,T_outc,T_inc)
    eff_new = hxf.effectiveness(heat_transfer_new,Cc,Ch,T_inc,T_inh)
    per_e_Q = (heat_transfer - heat_transfer_new)/heat_transfer
    per_e_eff = (eff - eff_new)/eff
    eff = eff_new
    heat_transfer = heat_transfer_new
    Q_counter += 1

    print(heat_transfer)


    if Q_counter > 100:
        print('exceeded max iterations for Q')
        break

    #now loop over entire thing again using these 4 values to get new property values to make answer more accurate. when this converges, can find effectiveness and Q
    #need to use lmtd and e-ntu approaches


def NTU():
    #ntu approach
    cmin = min(Cc,Ch)
    cmax = max(Cc,Ch)
    qmax = cmin * (T_inh - T_outc) #maximum possible heat tranfer
    Cr = cmin/cmax #ratio of specific heats 
    U = hxf.U_inside(hi,ho,hx.tube.d_inner,hx.tube.d_outer,hx.tube_length)
    A = 0.5 * np.pi*(hx.do * hx.L + hx.di * hx.L) * hx.tube_number #currently uses average of inside and outside area, check
    NTU = (U * A)/cmin
    if  hx.co_counter == 'counter':
        e = (1 - np.exp(-NTU * (1 + Cr)))/(1 + Cr) #equations from wiki, check
    if hx.co_counter == 'co':
        e = (1 - np.exp(-NTU * (1 - Cr)))/(1 - Cr * np.exp(-NTU * (1 - Cr)))
    else:
        print('Error, heat exchanger must be counter or co flow') 
    #may need something about mixed flow here later
    q = qmax * e
    return q
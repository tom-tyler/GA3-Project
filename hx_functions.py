from re import M
from tkinter import N
from scipy.optimize import fsolve
import numpy as np
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
from hx_classes import HX, water, pipe
import pandas as pd
from numpy import pi
from math import sqrt
import openpyxl

run1 = 1

#HYDRAULIC DESIGN

def hydraulic_design(m_c,m_h,h_w,c_w,hx,K_hot = 1.8,K_cold = 1):
    
    m_counter = 0
    dP_e_h = 1
    dP_e_c = 1

    #all class variables calculated here so that only need to do them once
    sigma = hx.sigma
    sigma_nozzle = hx.sigma_nozzle
    tube_number = hx.tube_number
    tube_passes = hx.tube_passes
    tube_c_area = hx.tube.c_area
    nozzle_c_area = hx.nozzle_c_area
    A_shell = hx.A_shell
    m_increment  = hx.m_increment
    pump_year = hx.pump_year
    accuracy = hx.accuracy

    for mdp_counter in range(100):
        if (abs(dP_e_h) < accuracy) and (abs(dP_e_c) < accuracy): 
            break
        #heat capacities
        Cc = m_c*c_w.cp
        Ch = m_h*h_w.cp

        #mass flow in one tube
        m_tube = m_h/(tube_number/tube_passes)

        #various velocities needed
        V_tube = V(m_tube,h_w,tube_c_area)
        V_nozzle_h = V(m_h,h_w,nozzle_c_area)
        V_shell = V(m_c,c_w,A_shell)
        V_nozzle_c = V(m_c,c_w,nozzle_c_area)

        
        #pressure drop in a single tube
        dP_tube = dP_tube_drop(hx,h_w,V_tube) * tube_passes

        #pressure drop due to separation leaving tubes
        dP_in_plus_out = dP_inout(h_w,V_tube,sigma) * tube_passes

        #pressure drop due to separation leaving nozzle
        dP_in_plus_out_nozzle = dP_inout(h_w,V_nozzle_h,sigma_nozzle)

        #head loss in nozzle
        dP_nozzles_h = 2 * dP_nozzle(V_nozzle_h,h_w)

        #overall pressure drop
        dP_tube_ovr = K_hot*(dP_tube + dP_in_plus_out + dP_in_plus_out_nozzle)  + dP_nozzles_h

        dP_new_h = mdot_dP(m_h,'h',h_w,pump_year)

        dP_e_h = (dP_new_h - dP_tube_ovr)/dP_new_h

        if dP_tube_ovr < dP_new_h:
            if abs(dP_e_h) < 1:
                m_h += abs(dP_e_h)*hx.m_increment/10
            else:
                m_h += hx.m_increment
        else:
            if abs(dP_e_h) < 1:
                m_h -= abs(dP_e_h)*hx.m_increment/10
            else:
                m_h -= hx.m_increment

        #cold side
        dP_shell = dP_shell_drop(c_w, m_c, hx, K_cold)
        dP_nozzles_c = 2 * dP_nozzle(V_nozzle_c,c_w)

        dP_shell_ovr = dP_shell + dP_nozzles_c

        dP_new_c = mdot_dP(m_c,'c',c_w,pump_year)

        dP_e_c = (dP_new_c - dP_shell_ovr)/dP_new_c

        if dP_shell_ovr < dP_new_c:
            if abs(dP_e_c) < 1:
                m_c += abs(dP_e_c)*hx.m_increment/10
            else:
                m_c += hx.m_increment
        else:
            if abs(dP_e_c) < 1:
                m_c -= abs(dP_e_c)*hx.m_increment/10
            else:
                m_c -= hx.m_increment

        m_counter += 1

        #print('m_h: ',m_h, '    m_c: ',m_c, '      dP_e_c: ',dP_e_c, '      dP_e_c: ',dP_e_c)
        if mdp_counter == 99:
            print('exceeded max iterations for m,dP')
            

    hydraulic = {'m_h':m_h,'m_c':m_c,'V_tube':V_tube,'V_shell':V_shell,'dP_hot':dP_tube_ovr,'dP_cold':dP_shell_ovr,'Cc':Cc,'Ch':Ch}
    return hydraulic



def dP_tube_drop(hx,liquid,V):

    di = hx.tube.d_inner
    L = hx.tube_length
    G = V*liquid.rho
    s = liquid.rho/1000

    f = friction_factor(Re(V,di,liquid))

    dP = f*(L/di)*0.5*liquid.rho*V**2
    #dP = (f*L*G**2)/(2000*di*s)    #same thing (same value)
    return dP



def Re(V,d,liquid):
    Re = (V*d*liquid.rho)/liquid.mu
    return Re



def dP_shell_drop(liquid, m_c, hx, K_cold):
    
    nb = hx.baffle_number

    dPi = dP_ideal(m_c, liquid, hx)
   
    dPw = dPw_ideal(m_c, liquid, hx)
  
    dPs = ((nb - 1)*dPi*hx.Rb + nb*dPw)*hx.Rl*K_cold * hx.shell_passes
    dPe = 2*dPi*(1 + hx.Ncw/hx.Nc)*hx.Rb*hx.Rs 
    dP = dPe + dPs
    return dP



def dPw_ideal(m_c, liquid, hx):
    dpw_ideal = ((2 + 0.6 * hx.Ncw) * m_c**2) / (2 * 1 * liquid.rho * hx.Sm * hx.Sw)
    return dpw_ideal



def dP_ideal(m_c, liquid, hx):
    #ideal pressure drop
    G = m_c/hx.Sm
    b = hx.b3 / (1 + 0.14 * ((hx.tube.d_outer * G)/liquid.mu)**hx.b4)
    phi = 1
    f_ideal = hx.b1 * (1.33/(hx.pitch/hx.tube.d_outer))**b * ((hx.tube.d_outer * G)/liquid.mu)**hx.b2
    dP_ideal = (2 * f_ideal * hx.Nc * G**2)/(1 * liquid.rho * phi)
    return dP_ideal



def friction_factor(Re, eD=0):
    A = eD/3.7+(6.7/Re)**0.9
    fo = 1/(-2*np.log10(eD/3.7-5.02/Re*np.log10(A)))**2

    if eD:
        f = fsolve(lambda x: 1/x**0.5+2.0*np.log10(eD/3.7+2.51/Re/x**0.5), fo)
    else:
        f = fsolve(lambda x: 1/x**0.5-2.0*np.log10(Re*x**0.5)+0.8, fo)
    return f[0]



def V(m_dot,liquid,area):
    rho = liquid.rho
    V = m_dot/(rho*area)
    return V



def h_inner(m_h, liquid, hx):
    m_tube = m_h/(hx.tube_number/hx.tube_passes)
    V_tube = V(m_tube,liquid,hx.tube.c_area)

    di = hx.tube.d_inner

    R = Re(V_tube,di,liquid)

    if R > 1e4:
        Nu = 0.023*R**0.8*liquid.Pr**0.3
    else:
        Nu = 0.116*(R**(2/3) - 125)*liquid.Pr**(1/3)*(1 + (di/hx.tube_length)**(2/3))

    hi = Nu*liquid.k/di

    return hi



def LMTD(T_1in,T_2in,T_1out,T_2out):
    lmtd = ((T_2in - T_1out) - (T_2out - T_1in)) / np.log((T_2in - T_1out) / (T_2out - T_1in))
    return lmtd



def U_inside(hi,ho,hx,k_copper = 398):

    do = hx.tube.d_outer
    di = hx.tube.d_inner
    L = hx.tube_length*hx.tube_passes

    ro = do/2
    ri = di/2
    Ao = np.pi*do*L
    Ai = np.pi*di*L

    inv_U = 1/hi + (Ai*np.log(ro/ri))/(2*np.pi*k_copper*L) + Ai/(Ao*ho)
    U = 1/inv_U
    return U



def KcKe(sigma):

    #values of exit and entrance pressure loss coefficient for turbulent flow, Re = 10 000
    Kc =  0.51 - (0.41/1) * sigma

    Ke_array = np.array([[0,1],
        [0.1,0.8],
        [0.2,0.62],
        [0.3,0.46],
        [0.4,0.32],
        [0.5,0.2],
        [0.6,0.1],
        [0.7,0.03],
        [0.8,-0.03],
        [0.9,-0.07],
        [1.0,-0.1]]) #sigma,Ke

    Ke_f = interp1d(Ke_array[:,0],Ke_array[:,1])
    
    Ke = Ke_f(sigma)

    return float(Kc),float(Ke)



def dP_inout(liquid,V,sigma):

    Kc,Ke = KcKe(sigma)

    dP = 0.5*liquid.rho*(V**2)*(Kc + Ke)

    return dP



def dP_nozzle(V,liquid):
    G = V*liquid.rho
    s = liquid.rho/1000
    

    dP = 7.5e-4*(G**2)/s
    return dP



def h_outer(m_c, liquid, hx): # bell delaware method

    #ideal calculation
    G = m_c / hx.Sm
    Re_S = (hx.tube.d_outer * G)/liquid.mu
    a = hx.a3 / (1 + 0.14 * (Re_S**hx.a4))
    j = hx.a1 * (1.33 /(hx.pitch/hx.tube.d_outer))**a * (Re_S**hx.a2) #dimensionless quanitity
    phi = 1 #viscosity correction factor
    ho_ideal = j * liquid.cp * G * liquid.Pr**(-2/3)*phi

    ho = ho_ideal * hx.Jc * hx.Jl * hx.Jb * hx.Jr * hx.Js
    return ho



def thermal_design(m_h,m_c,h_w,c_w,hx,T_inh,T_inc,T_outh,T_outc):
    
    #heat capacities
    Cc = m_c*c_w.cp
    Ch = m_h*h_w.cp

    hi = h_inner(m_h, h_w, hx)
    ho = h_outer(m_c, c_w, hx)
    
    U = U_inside(hi,ho,hx)
    A_con = hx.convection_area

    cmin = min(Cc,Ch)
    cmax = max(Cc,Ch)
    qmax = cmin * (T_inh - T_inc) #maximum possible heat tranfer
    Cr = cmin/cmax #ratio of specific heats 
    NTU = (U * A_con)/cmin
    c_root = (1 + Cr**2)**0.5
    NTU1 = NTU/hx.shell_passes
    e1 = fsolve(lambda e1: NTU1 + (np.log(((2/e1 - (1 + Cr))/c_root - 1)/((2/e1 - (1 + Cr))/c_root + 1)))/c_root, 0.1)[0]

    #e1 = 2 / (1 + Cr + c_root * ((1 + np.exp(-NTU*c_root))/(1 - np.exp(-NTU*c_root))))
    ez = ((1 - e1*Cr)/(1 - e1))**hx.shell_passes
    e = ((ez) - 1) / ((ez) - Cr)

    q_ntu = qmax * e

    T_outc = q_ntu/Cc + T_inc
    T_outh = T_inh - q_ntu/Ch


    #now solve thermal design equations by iteration to get T_outh and T_outc. also find P_outh and P_outc

    # T_outc_new,T_outh_new = 1,1
    # T_counter = 0
    # rel_e_h,rel_e_c = 1,1

    # while (abs(rel_e_c) > hx.accuracy) or (abs(rel_e_h) > hx.accuracy):

    #     F = correction_factor(T_inc,T_inh,T_outc,T_outh,hx)

    #     #might need to think about co/counterflow here for when shell passes > 1

    #     T_outc_new = fsolve(lambda T_outc: (T_inc - T_outc) + (1/Cc)*A_con*U*correction_factor(T_inc,T_inh,T_outc,T_outh,hx)*LMTD(T_inc,T_inh,T_outc,T_outh), T_outc)[0]
    
    #     T_outh_new = fsolve(lambda T_outh: (T_inh - T_outh) - (1/Ch)*A_con*U*correction_factor(T_inc,T_inh,T_outc,T_outh,hx)*LMTD(T_inc,T_inh,T_outc,T_outh), T_outh)[0]
        
    #     rel_e_c = (T_outc_new - T_outc)/T_outc
    #     rel_e_h = (T_outh_new - T_outh)/T_outh

    #     if T_outc > T_outc_new:
    #         T_outc = T_outc_new - hx.accuracy
    #     else:
    #         T_outc = T_outc_new + hx.accuracy

    #     if T_outh > T_outh_new:
    #         T_outh = T_outh_new - hx.accuracy
    #     else:
    #         T_outh = T_outh_new + hx.accuracy
    #     #T_outc = (T_outc_new + T_outc)/2
    #     #T_outh = (T_outh_new + T_outh)/2
    #     T_counter += 1

    #     if T_counter > 100:
    #         print('exceeded max iterations for T')
    #         break

    thermal = {'T_outh':T_outh,'T_outc':T_outc,'q_ntu':q_ntu,'eff_ntu':e,'U':U}
    return thermal



def Q_c(Cc,T_outc,T_inc):
    Q = Cc*(T_outc - T_inc)
    return Q



def Q_h(Ch,T_outh,T_inh):
    Q = Ch*(T_inh - T_outh)
    return Q



def effectiveness(Q,Cc,Ch,T_inc,T_inh):

    if Cc > Ch:
        e = Q / (Ch*(T_inh - T_inc))
    else:
        e = Q / (Cc*(T_inh - T_inc))
    
    return e



def mdot_dP(m_dot,side,liquid,year):

    if year == 2022:
        if side == 'h':
            mdot_dP_array = np.array([[0.4583,   0.1333e5],
                                    [0.4236,   0.1756e5],
                                    [0.4010,   0.2024e5],
                                    [0.3611,   0.2577e5],
                                    [0.3125,   0.3171e5],
                                    [0.2639,   0.3633e5],
                                    [0.2222,   0.4233e5],
                                    [0.1597,   0.4784e5],
                                    [0.1181,   0.5330e5],
                                    [0.0694,   0.5715e5]]) #[mass flow rate,pressure difference] for hot side
        elif side == 'c':
            mdot_dP_array = np.array([[0.5833,   0.1113e5],
                                    [0.5083,   0.2157e5],
                                    [0.4750,   0.2538e5],
                                    [0.4250,   0.3168e5],
                                    [0.3792,   0.3613e5],
                                    [0.3417,   0.4031e5],
                                    [0.2958,   0.4511e5],
                                    [0.2583,   0.4846e5],
                                    [0.2125,   0.5181e5],
                                    [0.1708,   0.5573e5]]) #[mass flow rate, pressure difference] for cold side

        else:
            print('please input side correctly')
    
    elif year == 2020:
        if side == 'h':
            mdot_dP_array = np.array([[0.5382,   0.1101e5],
                                      [0.5278,   0.1315e5],
                                      [0.4931,   0.1800e5],
                                      [0.4549,   0.2185e5],
                                      [0.4201,   0.2537e5],
                                      [0.3854,   0.2999e5],
                                      [0.3507,   0.3440e5],
                                      [0.3160,   0.3780e5],
                                      [0.2813,   0.4149e5],
                                      [0.2465,   0.4547e5],
                                      [0.2118,   0.5005e5],
                                      [0.1771,   0.5271e5],
                                      [0.1424,   0.5677e5],
                                      [0.1076,   0.5971e5],
                                      [0.0694,   0.6045e5]]) #[mass flow rate,pressure difference] for hot side
        elif side == 'c':
            mdot_dP_array = np.array([[0.6917,   0.1475e5],
                                      [0.6750,   0.1619e5],
                                      [0.6292,   0.2178e5],
                                      [0.5917,   0.2607e5],
                                      [0.5458,   0.3041e5],
                                      [0.5083,   0.3417e5],
                                      [0.4625,   0.3756e5],
                                      [0.4250,   0.4118e5],
                                      [0.3792,   0.4423e5],
                                      [0.3417,   0.4711e5],
                                      [0.2958,   0.5031e5],
                                      [0.2542,   0.5297e5],
                                      [0.2125,   0.5561e5],
                                      [0.1708,   0.5823e5]]) #[mass flow rate, pressure difference] for cold side

        else:
            print('please input side correctly')

    elif year == 2019:
        if side == 'h':
            mdot_dP_array = np.array([[0.5382,   0.1101e5],
                                      [0.5278,   0.1315e5],
                                      [0.4931,   0.1800e5],
                                      [0.4549,   0.2185e5],
                                      [0.4201,   0.2537e5],
                                      [0.3854,   0.2999e5],
                                      [0.3507,   0.3440e5],
                                      [0.3160,   0.3780e5],
                                      [0.2813,   0.4149e5],
                                      [0.2465,   0.4547e5],
                                      [0.2118,   0.5005e5],
                                      [0.1771,   0.5271e5],
                                      [0.1424,   0.5677e5],
                                      [0.1076,   0.5971e5],
                                      [0.0694,   0.6045e5]]) #[mass flow rate,pressure difference] for hot side
        elif side == 'c':
            mdot_dP_array = np.array([[0.6917,   0.1475e5],
                                      [0.6750,   0.1619e5],
                                      [0.6292,   0.2178e5],
                                      [0.5917,   0.2607e5],
                                      [0.5458,   0.3041e5],
                                      [0.5083,   0.3417e5],
                                      [0.4625,   0.3756e5],
                                      [0.4250,   0.4118e5],
                                      [0.3792,   0.4423e5],
                                      [0.3417,   0.4711e5],
                                      [0.2958,   0.5031e5],
                                      [0.2542,   0.5297e5],
                                      [0.2125,   0.5561e5],
                                      [0.1708,   0.5823e5]]) #[mass flow rate, pressure difference] for cold side

        else:
            print('please input side correctly')

    elif year == 2018:
        if side == 'h':
            mdot_dP_array = np.array([[0.4954,   0.0989e5],
                                      [0.4805,   0.1245e5],
                                      [0.4640,   0.1541e5],
                                      [0.4475,   0.1827e5],
                                      [0.4310,   0.2083e5],
                                      [0.4145,   0.2339e5],
                                      [0.3980,   0.2625e5],
                                      [0.3815,   0.2880e5],
                                      [0.3650,   0.3115e5],
                                      [0.3485,   0.3330e5],
                                      [0.3320,   0.3575e5],
                                      [0.3155,   0.3800e5],
                                      [0.2990,   0.4014e5],
                                      [0.2825,   0.4249e5],
                                      [0.2660,   0.4503e5],
                                      [0.2495,   0.4647e5],
                                      [0.2330,   0.4900e5],
                                      [0.2165,   0.5134e5],
                                      [0.2000,   0.5337e5],
                                      [0.1819,   0.5470e5],
                                      [0.1670,   0.5703e5],
                                      [0.1472,   0.5966e5],
                                      [0.1307,   0.6068e5],
                                      [0.1142,   0.6150e5],
                                      [0.1010,   0.6242e5],
                                      [0.0845,   0.6304e5],
                                      [0.0680,   0.6375e5],
                                      [0.0515,   0.6457e5]]) #[mass flow rate,pressure difference] for hot side
        elif side == 'c':
            mdot_dP_array = np.array([[0.4426,   0.1068e5],
                                      [0.4255,   0.1418e5],
                                      [0.4055,   0.1779e5],
                                      [0.3913,   0.2056e5],
                                      [0.3799,   0.2382e5],
                                      [0.3628,   0.2601e5],
                                      [0.3485,   0.2858e5],
                                      [0.3286,   0.3187e5],
                                      [0.3058,   0.3627e5],
                                      [0.2801,   0.4037e5],
                                      [0.2573,   0.4426e5],
                                      [0.2317,   0.4845e5],
                                      [0.2060,   0.5213e5],
                                      [0.1861,   0.5569e5],
                                      [0.1576,   0.6036e5],
                                      [0.1319,   0.6412e5],
                                      [0.1034,   0.6838e5],
                                      [0.0806,   0.7121e5],
                                      [0.0664,   0.7343e5],
                                      [0.0521,   0.7744e5]]) #[mass flow rate, pressure difference] for cold side

        else:
            print('please input side correctly')

    elif year == 2017:
        if side == 'h':
            mdot_dP_array = np.array([[0.4937,   0.0579e5],
                                      [0.4789,   0.0845e5],
                                      [0.4640,   0.1091e5],
                                      [0.4475,   0.1317e5],
                                      [0.4294,   0.1644e5],
                                      [0.4129,   0.1890e5],
                                      [0.3980,   0.2125e5],
                                      [0.3782,   0.2401e5],
                                      [0.3650,   0.2575e5],
                                      [0.3452,   0.2891e5],
                                      [0.3320,   0.3055e5],
                                      [0.3188,   0.3239e5],
                                      [0.3007,   0.3494e5],
                                      [0.2792,   0.3750e5],
                                      [0.2644,   0.3953e5],
                                      [0.2479,   0.4137e5],
                                      [0.2330,   0.4360e5],
                                      [0.2165,   0.4564e5],
                                      [0.2000,   0.4727e5],
                                      [0.1835,   0.4970e5],
                                      [0.1670,   0.5073e5],
                                      [0.1489,   0.5276e5],
                                      [0.1340,   0.5478e5],
                                      [0.1175,   0.5570e5],
                                      [0.0994,   0.5632e5],
                                      [0.0845,   0.5734e5],
                                      [0.0664,   0.5805e5]]) #[mass flow rate,pressure difference] for hot side
        elif side == 'c':
            mdot_dP_array = np.array([[0.4967,   0.0674e5],
                                      [0.4739,   0.1309e5],
                                      [0.4511,   0.2043e5],
                                      [0.4312,   0.2654e5],
                                      [0.4055,   0.3279e5],
                                      [0.3856,   0.3829e5],
                                      [0.3628,   0.4391e5],
                                      [0.3371,   0.4903e5],
                                      [0.3115,   0.5415e5],
                                      [0.2887,   0.5824e5],
                                      [0.2659,   0.6203e5],
                                      [0.2431,   0.6541e5],
                                      [0.2231,   0.6848e5],
                                      [0.2003,   0.7105e5],
                                      [0.1775,   0.7361e5],
                                      [0.1547,   0.7577e5],
                                      [0.1376,   0.7721e5],
                                      [0.1120,   0.7916e5],
                                      [0.0664,   0.8183e5]]) #[mass flow rate, pressure difference] for cold side

        else:
            print('please input side correctly')

    else:
        print('please input a year')

    # if dP_ovr > 0.65e5:
    #     dP_ovr = 0.65e5

    #mdot_from_dP = interp1d(mdot_dP_array[:,1],mdot_dP_array[:,0],fill_value='extrapolate',kind = 'linear')
    dP_from_mdot = interp1d(mdot_dP_array[:,0],mdot_dP_array[:,1],fill_value='extrapolate',kind = 'cubic')

    dP_new = dP_from_mdot(m_dot/(liquid.rho/1000))
    #m_dot = mdot_from_dP(dP_ovr)*(liquid.rho/1000)

    #rel_e_dP = (dP_new - dP_ovr)/dP_ovr
    #rel_e_mdot = (m_dot_new - m_dot)/m_dot

    #plt.plot(mdot_dP_array[:,0],mdot_dP_array[:,1])
    #plt.xlabel('m')
    #plt.ylabel('p')

    return dP_new



def correction_factor(T_inc,T_inh,T_outc,T_outh,hx):

    N = hx.shell_passes
    #N = number of shell passes. 2M = tube passes per shell

    P = (T_outh - T_inh)/(T_inc - T_inh)

    R = (T_inc - T_outc)/(T_outh - T_inh)

    S = ((R**2 + 1)**0.5)/(R - 1)

    W = ((1-P*R)/(1 - P))**(1/N)

    F = (S*np.log(W))/np.log((1 + W - S + S*W)/(1 + W + S - S*W))

    return F



def hx_design(hx,K_hot,K_cold):

    #initial guesses for mass flowrate:
    m_h = 0.45 #initial guess for hot mass flow rate
    m_c = 0.45 #initial guess for cold mass flow rate

    #initial guesses for outlet temperatures
    T_outh = hx.T_inh - 9.103
    T_outc = hx.T_inc + 10.645
    R = (hx.T_inc - T_outc)/(T_outh - hx.T_inh)

    if (abs(1-R) <0.1):
        T_outh = hx.T_inh - 6.2
        T_outc = hx.T_inc + 7.876


    #creating hot and cold water objects
    if (T_outh > 70) or (T_outh < 30):
        T_outh = 50
    if (T_outc > 40) or (T_outc < 10):
        T_outc = 20
    h_w = water(hx.T_inh,T_outh)
    c_w = water(hx.T_inc,T_outc)

    #HYDRAULIC DESIGN
    hydraulic = hydraulic_design(m_c,m_h,h_w,c_w,hx,K_hot,K_cold) #,invalid_hx_flag 


    m_h, m_c = hydraulic['m_h'], hydraulic['m_c']
    dP_hot, dP_cold = hydraulic['dP_hot'], hydraulic['dP_cold']

    #THERMAL DESIGN
    thermal = thermal_design(m_h,m_c,h_w,c_w,hx,hx.T_inh,hx.T_inc,T_outh,T_outc)
    T_outh, T_outc = thermal['T_outh'], thermal['T_outc']

    heat_transfer_ntu = thermal['q_ntu']
    eff_ntu = thermal['eff_ntu']
    U = thermal['U']

    hx_dict = {'Name':f'{hx.name}-p{(hx.pump_year)}',
            'T cold in (C)':hx.T_inc,
            'T cold out (C)':T_outc,
            'T hot in (C)':hx.T_inh,
            'T hot out (C)':T_outh,
            'mdot_cold (l/s)':m_c/(c_w.rho/1000),
            'mdot_hot (l/s)':m_h/(h_w.rho/1000),
            'dP_cold (bar)':dP_cold/1e5,
            'dP_hot (bar)':dP_hot/1e5,
            'Q_NTU (kW)':heat_transfer_ntu,
            'eff_NTU':eff_ntu,
            'mass (kg)':hx.total_mass()
            }                                                        
            
    return hx_dict


def heat_exchangers(heat_exchanger=None):

    hx_list = {}

    hx_list['hx_y2018_p2022'] = HX(tube_number = 13,
                        baffle_number = 14,
                        pitch = 12e-3,
                        tube_length = 362e-3,
                        plenum_length_1 = 50e-3,
                        plenum_length_2 = 50e-3,
                        baffle_gap = 14e-3,
                        baffle_type = 'across_c',
                        tube_layout='t',
                        shell_passes = 1,
                        tube_bundle_diameter= 56e-3,
                        tube_passes = 1,
                        design_year = 2018,
                        pump_year = 2022,
                        T_inh = 53.4,
                        T_inc = 19.7,
                        baffle_spacing_in = 39e-3,
                        baffle_spacing_out = 39e-3,
                        leakage = True,
                        name = 'JPL-2018'
                        )


    hx_list['hx_y2018_p2019'] = HX(tube_number = 13,
                        baffle_number = 14,
                        pitch = 12e-3,
                        tube_length = 362e-3,
                        plenum_length_1 = 50e-3,
                        plenum_length_2 = 50e-3,
                        baffle_gap = 14e-3,
                        baffle_type = 'across_c',
                        tube_layout='t',
                        shell_passes = 1,
                        tube_bundle_diameter= 56e-3,
                        tube_passes = 1,
                        baffle_spacing_in = 39e-3,
                        baffle_spacing_out = 39e-3,
                        design_year = 2018,
                        pump_year = 2019,
                        T_inh = 53.2,
                        T_inc = 13.8,
                        leakage = True,
                        name = 'JPL-2018')

    hx_list['hx_y2017B_p2022'] = HX(tube_number = 14,
                        baffle_number = 12,
                        pitch = 15e-3,
                        tube_length = 225e-3,
                        plenum_length_1 = 41e-3,
                        plenum_length_2 = 23e-3,
                        baffle_gap = 25.6e-3,
                        baffle_type = 'across_c',
                        tube_layout='t',
                        shell_passes = 2,
                        tube_bundle_diameter= 60e-3,
                        tube_passes = 2,
                        baffle_spacing_in = 31e-3,
                        baffle_spacing_out = 33e-3,
                        design_year = 2017,
                        pump_year = 2022,
                        T_inh = 48.5,
                        T_inc = 22.3,
                        leakage = True,
                        name = '2017-B')

    hx_list['hx_y2017B_p2017'] = HX(tube_number = 14,
                        baffle_number = 12,
                        pitch = 15e-3,
                        tube_length = 225e-3,
                        plenum_length_1 = 41e-3,
                        plenum_length_2 = 23e-3,
                        baffle_gap = 25.6e-3,
                        baffle_type = 'across_c',
                        tube_layout='t',
                        shell_passes = 2,
                        tube_bundle_diameter= 60e-3,
                        tube_passes = 2,
                        baffle_spacing_in = 31e-3,
                        baffle_spacing_out = 33e-3,
                        design_year = 2017,
                        pump_year = 2017,
                        T_inh = 52.4,
                        T_inc = 21.7,
                        leakage = True,
                        name = '2017-B')

    hx_list['hx_y2017C_p2017'] = HX(tube_number = 14,
                        baffle_number = 8,
                        pitch = 10e-3,
                        tube_length = 192e-3,
                        plenum_length_1 = 41e-3,
                        plenum_length_2 = 24e-3,
                        baffle_gap = 15.75e-3,
                        baffle_type = 'across_c',
                        tube_layout='t',
                        shell_passes = 1,
                        tube_bundle_diameter= 52e-3,
                        tube_passes = 2,
                        baffle_spacing_in = 31.25e-3,
                        baffle_spacing_out = 31.25e-3,
                        design_year = 2017,
                        pump_year = 2017,
                        T_inh = 51.5,
                        T_inc = 22.2,
                        leakage = True,
                        name = '2017-C')

    if heat_exchanger == None:
        return hx_list
    else:
        hx_singular = {}
        hx_singular[heat_exchanger] = hx_list[heat_exchanger]
        return hx_singular



def hx_design_basic(hx):
    #initial guesses for mass flowrate:
    m_h = 0.45 #initial guess for hot mass flow rate
    m_c = 0.45 #initial guess for cold mass flow rate

    #initial guesses for outlet temperatures
    T_outh = hx.T_inh - 9.103
    T_outc = hx.T_inc + 10.645
    R = (hx.T_inc - T_outc)/(T_outh - hx.T_inh)

    if (abs(1-R) <0.1):
        T_outh = hx.T_inh - 6.2
        T_outc = hx.T_inc + 7.876

    #initial heat transfer parameters
    heat_transfer,eff = 1,1
    Q_counter = 0
    rel_e_h1,rel_e_c1 = 1,1
    q_acc = 0.05



def hx_moodle_data(K_hot = 1.8, K_cold = 1):


    hx_list = heat_exchangers()
    hx_data = pd.DataFrame()

    for hxi in hx_list:
        hx = hx_list[hxi]

        performance = hx_design(hx,K_hot,K_cold) #,invalid_hx_flag 


        hx_data = hx_data.append(performance, ignore_index = True)

    #order columns nicely
    hx_data = hx_data[['Name',
                'T cold in (C)',
                'T cold out (C)',
                'T hot in (C)',
                'T hot out (C)',
                'mdot_cold (l/s)',
                'mdot_hot (l/s)',
                'dP_cold (bar)',
                'dP_hot (bar)',
                'Q_NTU (kW)',
                'eff_NTU',
                'mass (kg)'
                    ]]
    hx_data = hx_data.sort_values(by="Q_NTU (kW)", ascending=False).head()
    with pd.option_context('display.max_rows', None, 'display.max_columns', None,"display.precision", 3):  # more options can be specified also
        print(hx_data)



def brute_opt(n = 10,K_hot = 1.8,K_cold = 1):
    #brute force optimisation with a few checks to eliminate cases as early as possible

    baffle_type = 'across_c' #assume this is best
    tube_layout = 't'        #assume this is best

    hx_designs = {}
    hx_data = pd.DataFrame()
    
    design_no = 0

    tp_array = np.array([4])
    sp_array = np.array([1])
    pl1_array = np.array([41e-3]) #np.linspace(41e-3,100e-3,n)
    bsi_array = np.array([41e-3]) #np.linspace(41e-3,50e-3,n)
    bn_array = np.array(range(0,10 + 1))   #
    bg_array = np.linspace(10e-3,30e-3,int(n/2))
    p_array = np.linspace(10e-3,14e-3,int(n/2))

    packing_density = np.pi/np.sqrt(12)
    r = 64e-3/2
    segment_area = r**2 * np.arccos(0.8660) - (0.8660*r)*(2*r*0.13397*r - (0.13397*r)**2)**0.5
    max_tube_area = (np.pi*(r**2) - 6*segment_area)*packing_density

    shortest_tube = 100e-3

    tube = pipe(6e-3,8e-3,0.20,3.5,shortest_tube)

    for pitch in p_array:
        print('pitch: ',pitch)

        eff_tube_area = np.pi*(pitch/2)**2
        max_no_tubes_from_area = int(max_tube_area/eff_tube_area)
        max_no_tube_from_mass = int(1.1/tube.mass)
        max_no_tubes = min(max_no_tubes_from_area,max_no_tube_from_mass)
        min_no_tubes = int(max_no_tubes*0.4)
        tn_array = np.array(range(min_no_tubes,max_no_tubes + 1))

        for tube_passes in tp_array:
            #print('tube_passes: ',tube_passes)

            # if tube_passes%2 == 0:
            #     l_min = 10e-3
            # else:
            #     l_min = 41e-3
            pl2_array = np.array([41e-3]) #np.linspace(l_min,100e-3,n)
        
            for shell_passes in sp_array:
                #print('shell_passes: ',shell_passes)

                # if shell_passes%2 == 0:
                #     bso_min = 10e-3
                # else:
                #     bso_min = 41e-3
                bso_array = np.array([41e-3]) #np.linspace(bso_min,50e-3,n)
                for tube_number in tn_array:
                    #print('tube_number',tube_number)

                    for plenum_length_1 in pl1_array:                        

                        for plenum_length_2 in pl2_array:
                    
                            max_tube_length = (0.35 - plenum_length_1 - plenum_length_2)
                            tl_array = np.linspace(shortest_tube,max_tube_length,n)
                            n_total = tp_array.size*sp_array.size*tl_array.size*pl1_array.size*bsi_array.size*bn_array.size*bg_array.size*p_array.size*bso_array.size*tn_array.size*pl2_array.size

                            for tube_length in tl_array:

                                tube.l = tube_length

                                if (tube_number * tube_length <= 3.5) and (tube_number * tube.mass <= 1.1):               
                                
                                    for baffle_spacing_in in bsi_array:

                                        for baffle_spacing_out in bso_array:

                                            for baffle_number in bn_array:
                                                #print('baffle_number',baffle_number)

                                                for baffle_gap in bg_array:    
                                                                
                                                    design_no += 1
                                                    #print(f'design number = {design_no}, designs = {n_total}')
                                                    heat_exchanger = HX(tube_number = tube_number,
                                                                        baffle_number = baffle_number,
                                                                        pitch = pitch,
                                                                        tube_length = tube_length,
                                                                        plenum_length_1 = plenum_length_1,
                                                                        plenum_length_2 = plenum_length_2,
                                                                        baffle_gap = baffle_gap,
                                                                        baffle_type = baffle_type,
                                                                        tube_layout = tube_layout,
                                                                        shell_passes = shell_passes,
                                                                        tube_bundle_diameter  = (64e-3- (pitch-tube.d_outer)),
                                                                        tube_passes = tube_passes,
                                                                        baffle_spacing_in = baffle_spacing_in,
                                                                        baffle_spacing_out = baffle_spacing_out,
                                                                        design_year = 2022,
                                                                        pump_year = 2022,
                                                                        T_inh = 53.4,
                                                                        T_inc = 19.2,
                                                                        leakage = True,
                                                                        name = f'design number = {design_no}',
                                                                        co_counter='counter',
                                                                        approximate_glue_mass=0
                                                                        )
                                                    
                                                    if heat_exchanger.total_mass() <= 1.1:
                                                        hx_designs[f'design {design_no}'] = heat_exchanger
                                                        performance = hx_design(heat_exchanger,K_hot,K_cold) #,invalid_hx_flag
                                                        design = vars(heat_exchanger)
                                                        performance.update(design)
                                                        #if invalid_hx_flag == False:
                                                        hx_data = hx_data.append(performance, design) 
                                                            

    #order columns nicely
    hx_data = hx_data.sort_values(by="Q_NTU (kW)", ascending=False).head(10000)[['Name',
                                                                                'Q_NTU (kW)',
                                                                                'eff_NTU',
                                                                                'mass (kg)',
                                                                                'tube_number',
                                                                                'baffle_number',
                                                                                'pitch',
                                                                                'tube_length',
                                                                                'plenum_length_1',
                                                                                'plenum_length_2',
                                                                                'baffle_gap',
                                                                                'baffle_type',
                                                                                'tube_layout',
                                                                                'shell_passes',
                                                                                'tube_bundle_diameter',
                                                                                'tube_passes',
                                                                                'baffle_spacing_in',
                                                                                'baffle_spacing_out'
                                                                                ]]
    hx_data = hx_data
    with pd.option_context('display.max_rows', None, 'display.max_columns', None,"display.precision", 3):  # more options can be specified also
        print(hx_data)

    hx_data.to_excel(f"hx_data_{run1}.xlsx", sheet_name="heat_exchanger_data", index=False)

def brute_opt_2():
    #brute optimisation but applying some common sense to reduce time
    baffle_type = 'across_c'
    tube_layout = 't'

    hx_designs = {}
    hx_data = []
    K_hot = 1.8
    K_cold = 1
    
    design_no = 0

    pitch = 12e-3
    crossflow_rows = 4.5
    tube_bundle_diameter = (crossflow_rows * pitch + 8e-3)/1000
    plenum_length_1 = 41
    plenum_length_2 = 41

    for tube_passes in range(1,4):
        
        for tube_number in range(10,15):
            
            for tube_length in range(100,250,30):
                if tube_number * tube_length <= 3500:
                    for shell_passes in range(1,2):
                        if shell_passes%2 == 0:
                            bso_min = 10
                        else:
                            bso_min = 41
                        baffle_spacing_in = 41
                        for baffle_spacing_out in range(bso_min, 90, 40):
                            for baffle_number in range(10,14):
                                for baffle_gap in range(10,60,10):
                                    design_no += 1
                                    heat_exchanger = HX(tube_number = tube_number,
                                                        baffle_number = baffle_number,
                                                        pitch = 12/1000,
                                                        tube_length = tube_length/1000,
                                                        plenum_length_1 = plenum_length_1,
                                                        plenum_length_2 = plenum_length_2,
                                                        baffle_gap = baffle_gap/1000,
                                                        baffle_type = baffle_type,
                                                        tube_layout = tube_layout,
                                                        shell_passes = shell_passes,
                                                        tube_bundle_diameter  = tube_bundle_diameter,
                                                        tube_passes = tube_passes,
                                                        baffle_spacing_in = baffle_spacing_in/1000,
                                                        baffle_spacing_out = baffle_spacing_out/1000,
                                                        design_year = 2022,
                                                        pump_year = 2022,
                                                        T_inh = 53.4,
                                                        T_inc = 19.2,
                                                        leakage = True,
                                                        name = None,
                                                        co_counter='counter',
                                                        approximate_glue_mass=0
                                                        )
                                    
                                    if heat_exchanger.total_mass() <= 1.1:
                                        print('pass')
                                        hx_designs[f'design {design_no}'] = heat_exchanger
                                        performance = hx_design_basic(heat_exchanger)
                                        
                                        hx_data = hx_data.append(performance, vars(heat_exchanger)[0:13])
    
    hx_data.sort(reverse = True)
    return hx_data[0:3]

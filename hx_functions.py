from re import M
from tkinter import N
from scipy.optimize import fsolve
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from hx_classes import HX, water, pipe
import pandas as pd
from datetime import date



# PRESSURE DROP FUNCTIONS
#region

def dP_tube_drop(hx,liquid,V):
    '''Function returns the pressure drop inside a tube, given the tube velocity, liquid, and the heat exhanger'''
    di = hx.tube.d_inner
    L = hx.tube_length
    R = Re(V,di,liquid)
    f = friction_factor(R)

    #dP = f*(L/di)*0.5*liquid.rho*V**2
    dP = (4*f*L*hx.tube_passes/di + 4*hx.tube_passes)*liquid.rho*V**2/2
    return dP



def dP_shell_drop(liquid, m_c, hx,k4=1,k5=1):
    
    nb = hx.baffle_number
    G = m_c/hx.Sm
    b = hx.b3 / (1 + 0.14 * ((hx.tube.d_outer * G)/liquid.mu)**hx.b4)
   
    f_ideal = hx.b1 * (1.33/(hx.pitch/hx.tube.d_outer))**b * ((hx.tube.d_outer * G)/liquid.mu)**hx.b2

    dPi = k4 * 2*f_ideal*(G**2/liquid.rho)*hx.Nc
   
    dPw = k5 * (nb*hx.Rl*(2 + 0.6 * hx.Ncw)*G**2)/(2*liquid.rho)
  
    dPc = (nb - 1)*dPi*hx.Rl*hx.Rb

    dPe = dPi*(1 + hx.Ncw/hx.Nc)*hx.Rb*hx.Rs 

    dP_single_shell_pass = dPw + dPc + dPe
    dP = dP_single_shell_pass * hx.shell_passes
    return dP



def dP_inout(liquid,V,sigma):

    Kc,Ke = KcKe(sigma)

    dP = 0.5*liquid.rho*(V**2)*(Kc + Ke)

    return dP



def dP_nozzle(V,liquid):

    dP = 0.5*liquid.rho*V**2
    return dP

#endregion

# OFTEN USED VARIABLES
#region

def Re(V,d,liquid):
    R = (V*d*liquid.rho)/liquid.mu
    return R



def V(m_dot,liquid,area):
    rho = liquid.rho
    V = m_dot/(rho*area)
    return V

#endregion

# CORRELATIONS
#region
def friction_factor(R, eD=0):
    A = eD/3.7+(6.7/R)**0.9
    fo = 1/(-2*np.log10(eD/3.7-5.02/R*np.log10(A)))**2

    if eD:
        f = fsolve(lambda x: 1/x**0.5+2.0*np.log10(eD/3.7+2.51/R/x**0.5), fo)
    else:
        f = fsolve(lambda x: 1/x**0.5-2.0*np.log10(R*x**0.5)+0.8, fo)
    return f[0]


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
    #print(mdot_dP_array[-1,0])
    if m_dot < mdot_dP_array[-1,0]:
        m_dot = mdot_dP_array[-1,0] + 0.01

    #mdot_from_dP = interp1d(mdot_dP_array[:,1],mdot_dP_array[:,0],fill_value='extrapolate',kind = 'linear')
    dP_from_mdot = interp1d(mdot_dP_array[:,0],mdot_dP_array[:,1],fill_value='extrapolate',kind = 'cubic')

    dP_new = float(dP_from_mdot(m_dot/(liquid.rho/1000)))
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

#endregion

# HEAT TRANSFER VARIABLES
#region
def h_inner(m_h, liquid, hx, z1=1):
    m_tube = m_h/(hx.tube_number/hx.tube_passes)
    V_tube = V(m_tube,liquid,hx.tube.c_area)

    di = hx.tube.d_inner

    R = Re(V_tube,di,liquid)
    fr = z1*friction_factor(R)

    Nu = (fr/2)*R*liquid.Pr/(1.07 + 12.7*np.sqrt(fr/2)*(liquid.Pr**(2/3) - 1))

    hi = Nu*liquid.k/di

    return hi



def h_outer(m_c, liquid, hx, z2=1): # bell delaware method

    #ideal calculation
    G = m_c / hx.Sm
    Re_S = (hx.tube.d_outer * G)/liquid.mu

    a = hx.a3 / (1 + 0.14 * (Re_S**hx.a4))
    j = hx.a1 * (1.33 /(hx.pitch/hx.tube.d_outer))**a * (Re_S**hx.a2) #dimensionless quanitity

    ho_ideal = j * liquid.cp * G * liquid.Pr**(-2/3)

    ho = ho_ideal * hx.Jc * hx.Jl * hx.Jb * hx.Jr * hx.Js * z2
    return ho



def U_inside(hi,ho,hx,k_copper = 398):

    do = hx.tube.d_outer
    di = hx.tube.d_inner

    inv_U = 1/ho + do*np.log(do/di)/(2**k_copper) + do/(hi*di)
    U = 1/inv_U
    return U

#endregion

# DESIGN TOOLS
#region

# THERMAL DESIGN

def thermal_design(m_h,m_c,h_w,c_w,hx,z1=1,z2=1,z3=1):
    
    #m_c = z1*m_c
    #m_h = z2*m_h

    #heat capacities
    Cc = m_c*c_w.cp
    Ch = m_h*h_w.cp

    hi = h_inner(m_h, h_w, hx,z1=z1)
    ho = h_outer(m_c, c_w, hx,z2=z2)
    
    U = U_inside(hi,ho,hx)
    A_con = z3*hx.convection_area

    cmin = min(Cc,Ch)
    cmax = max(Cc,Ch)
    qmax = cmin * (hx.T_inh - hx.T_inc) #maximum possible heat tranfer
    Cr = cmin/cmax #ratio of specific heats 
    NTU = (U * A_con)/cmin
    c_root = (1 + Cr**2)**0.5

    e1 = 2 / (1 + Cr + c_root * ((1 + np.exp(-NTU*c_root))/(1 - np.exp(-NTU*c_root))))
    if hx.shell_passes == 1:
        e = e1
    else:
        ez = ((1 - e1*Cr)/(1 - e1))**(hx.shell_passes)
        e =((ez) - 1) / ((ez) - Cr)

    q_ntu = qmax * e
    q_corr = cmin * 40 * e

    T_outc = q_ntu/Cc + hx.T_inc
    T_outh = hx.T_inh - q_ntu/Ch

    thermal = {'T_outh':T_outh,'T_outc':T_outc,'q_ntu':q_ntu,'eff_ntu':e,'U':U,'q_corr':q_corr}
    return thermal

# HYDRAULIC DESIGN

def hydraulic_design(m_c,m_h,h_w,c_w,hx,k1=1,k2=1,k3=1,k4=1,k5=1,k6=1):
    
    mdp_counter = 0
    dP_e_h = 1
    dP_e_c = 1

    #all class variables calculated here so that only need to do them once
    sigma = hx.sigma
    tube_number = hx.tube_number
    tube_passes = hx.tube_passes
    tube_c_area = hx.tube.c_area
    nozzle_c_area = hx.nozzle_c_area
    A_shell = hx.A_shell
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

        #hot side loss
        #region
        #pressure drop in a single tube
        dP_tube = dP_tube_drop(hx,h_w,V_tube)

        #pressure drop due to separation leaving tubes
        dP_in_plus_out = dP_inout(h_w,V_tube,sigma) * tube_passes

        #head loss in nozzles
        dP_nozzles_h = 2 * dP_nozzle(V_nozzle_h,h_w)

        #overall pressure drop on hot side
        dP_tube_ovr = k1*dP_tube + k2*dP_in_plus_out + k3*dP_nozzles_h

        dP_new_h = mdot_dP(m_h,'h',h_w,pump_year)

        dP_e_h = (dP_new_h - dP_tube_ovr)/dP_new_h

        if dP_tube_ovr < dP_new_h:
            if abs(dP_e_h) < 1:
                m_h += abs(dP_e_h)*hx.m_increment*0.9
            else:
                m_h += hx.m_increment
        else:
            if abs(dP_e_h) < 1:
                m_h -= abs(dP_e_h)*hx.m_increment
            else:
                m_h -= hx.m_increment
        #endregion

        #cold side loss
        #region
        #shell pressure drop
        dP_shell = dP_shell_drop(c_w, m_c, hx,k4=k4,k5=k5)

        #head loss in nozzles
        dP_nozzles_c = 2 * dP_nozzle(V_nozzle_c,c_w)

        #overall pressure drop on cold side
        dP_shell_ovr = dP_shell + k6*dP_nozzles_c

        dP_new_c = mdot_dP(m_c,'c',c_w,pump_year)

        dP_e_c = (dP_new_c - dP_shell_ovr)/dP_new_c

        if dP_shell_ovr < dP_new_c:
            if abs(dP_e_c) < 1:
                m_c += abs(dP_e_c)*hx.m_increment*0.9
            else:
                m_c += hx.m_increment
        else:
            if abs(dP_e_c) < 1:
                m_c -= abs(dP_e_c)*hx.m_increment
            else:
                m_c -= hx.m_increment
        #endregion
        
        mdp_counter += 1

        if mdp_counter >= 10:
            pass
        elif mdp_counter >= 20: 
            hx.m_increment = 0.12
        elif mdp_counter >= 30:
            hx.m_increment = 0.9
        elif mdp_counter >= 99:
            print('exceeded max iterations for m,dP')
            
    hx.m_increment = 0.1
    hydraulic = {'m_h':m_h,'m_c':m_c,'V_tube':V_tube,'V_shell':V_shell,'dP_hot':dP_tube_ovr,'dP_cold':dP_shell_ovr,'Cc':Cc,'Ch':Ch}
    return hydraulic

# HX DESIGN: tool which outputs performance for a given heat exchanger in a dictionary

def hx_design(hx,k_array = np.ones(10)):

    #curve fitting constants
    k1 = k_array[0]
    k2 = k_array[1]
    k3 = k_array[2]
    k4 = k_array[3]
    k5 = k_array[4]
    k6 = k_array[5]
    z1 = k_array[6]
    z2 = k_array[7]
    z3 = k_array[8]
    mk = k_array[9]

    #initial guesses for mass flowrate:
    m_h = 0.45 #initial guess for hot mass flow rate
    m_c = 0.45 #initial guess for cold mass flow rate

    #initial guesses for outlet temperatures
    T_outh = hx.T_inh - 9.103
    T_outc = hx.T_inc + 10.645
    R = (hx.T_inc - T_outc)/(T_outh - hx.T_inh)

    if (abs(1-R) < 0.1):
        T_outh = hx.T_inh - 6.2
        T_outc = hx.T_inc + 7.876

    h_w = water(hx.T_inh,T_outh)
    c_w = water(hx.T_inc,T_outc)

    #HYDRAULIC DESIGN
    hydraulic = hydraulic_design(m_c,m_h,h_w,c_w,hx,k1=k1,k2=k2,k3=k3,k4=k4,k5=k5,k6=k6) #,invalid_hx_flag 

    m_h, m_c = hydraulic['m_h'], hydraulic['m_c']
    dP_hot, dP_cold = hydraulic['dP_hot'], hydraulic['dP_cold']
    Ch, Cc = hydraulic['Ch'], hydraulic['Cc']

    #THERMAL DESIGN
    thermal = thermal_design(m_h,m_c,h_w,c_w,hx,z1,z2,z3)
    T_outh, T_outc = thermal['T_outh'], thermal['T_outc']

    heat_transfer_ntu = thermal['q_ntu']
    eff_ntu = thermal['eff_ntu']
    ht_corr = thermal['q_corr']

    hx_dict = {'Name':f'{hx.name}-p{(hx.pump_year)}',
            'T cold in (C)':hx.T_inc,
            'T cold out (C)':T_outc,
            'T hot in (C)':hx.T_inh,
            'T hot out (C)':T_outh,
            'mdot_cold (l/s)':m_c/(c_w.rho/1000),
            'mdot_hot (l/s)':m_h/(h_w.rho/1000),
            'dP_cold (bar)':dP_cold/1e5,
            'dP_hot (bar)':dP_hot/1e5,
            'Q_NTU (kW)':heat_transfer_ntu/1000,
            'Q_NTU_corr (kW)':ht_corr/1000,
            'eff_NTU':eff_ntu,
            'mass (kg)':mk*hx.total_mass()
            }                                                        
            
    return hx_dict

#endregion

# CURVE FITTING
#region

def f_m_hot(i_list,k1,k2,k3):
    hx_dict = heat_exchanger_dict()
    k_array = np.array([k1,k2,k3,1,1,1,1,1,1,1])
    m_hot_list = []
    for i in i_list:
        hx = hx_dict[i]
        m_hot_list.append(hx_design(hx,k_array)['mdot_hot (l/s)'])

    return np.array(m_hot_list)

def f_m_cold(i_list,k4,k5,k6):
    hx_dict = heat_exchanger_dict()
    k_array = np.array([1,1,1,k4,k5,k6,1,1,1,1])
    m_cold_list = []
    for i in i_list:
        hx = hx_dict[i]
        m_cold_list.append(hx_design(hx,k_array)['mdot_cold (l/s)'])

    return np.array(m_cold_list)

def f_mass(i_list,mk):
    hx_dict = heat_exchanger_dict()
    k_array = np.array([1,1,1,1,1,1,1,1,1,mk])
    mass_predict_list = []
    for i in i_list:
        hx = hx_dict[i]
        mass_predict_list.append(hx_design(hx,k_array)['mass (kg)'])

    return np.array(mass_predict_list)

def fit_data(heat_exchanger=None):


    hx_dict = heat_exchanger_dict(heat_exchanger)
    print(hx_dict)
    hx_list2 = []
    m_list_hot = np.array([])
    m_list_cold = np.array([])
    q_list = np.array([])
    mass_list = np.array([])
    i_list = np.array(list(hx_dict.keys()))

    for hxi in hx_dict:
        hx = hx_dict[hxi]
        print(hx.name)

        experimental_data = hx.real_data
        m_list_hot = np.append(m_list_hot, np.array([experimental_data['mdot_hot (l/s)']]))
        m_list_cold = np.append(m_list_cold, np.array([experimental_data['mdot_cold (l/s)']]))
        q_list = np.append(q_list, np.array([experimental_data['Q_NTU_corr (kW)']]))
        mass_list = np.append(mass_list,np.array([experimental_data['mass (kg)']]))

    print('finding k1 k2 k3')
    k1k2k3, pcov1 = curve_fit(f_m_hot,i_list,m_list_hot,bounds=(0.75,2),p0 = np.array([1,1,1]))
    print('finding k4 k5 k6')
    k4k5k6, pcov2 = curve_fit(f_m_cold,i_list,m_list_cold,bounds=([0.5,0.5,0.75],[10,10,1.25]),p0 = np.array([1,1,1]))

    def f_thermal(i_list,z1,z2,z3):
        hx_dict = heat_exchanger_dict()
        k_array = np.array([k1k2k3[0],k1k2k3[1],k1k2k3[2],k4k5k6[0],k4k5k6[1],k4k5k6[2],z1,z2,z3,1])
        Q_predict_list = []
        for i in i_list:
            hx = hx_dict[i]
            Q_predict_list.append(hx_design(hx,k_array)['Q_NTU_corr (kW)'])

        return np.array(Q_predict_list)

    print('finding z1 z2 z3')
    z1z2z3, pcov3 = curve_fit(f_thermal,i_list,q_list,bounds=(0.5,1.5),p0 = np.array([1,1,1]))
    
    #z1z2z3 = np.array([1,1,1])
    
    print('findig mk')
    mk,pcov4 = curve_fit(f_mass,i_list,mass_list,bounds=(0.2,2),p0 = np.array(0.9))

    k_array = np.append(k1k2k3,k4k5k6)
    k_array = np.append(k_array,z1z2z3)
    k_array = np.append(k_array,mk)
    print('k1,k2,k3,k4,k5,k6,z1,z2,z3,mk = ',k_array)
    return k_array

#endregion

# DICTIONARY GENERATION
#region


def heat_exchanger_dict(heat_exchanger=None):

    hx_list = {}
    real_moodle_data = moodle_performance_dict()

    hx_list[1] = HX(tube_number = 13,
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
                        name = 'JPL-2018',
                        real_data = real_moodle_data['hx_y2018_p2022'])

    hx_list[2] = HX(tube_number = 14,
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
                        name = '2017-B',
                        real_data = real_moodle_data['hx_y2017B_p2022'])

    hx_list[3] = HX(tube_number = 13,
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
                        name = 'JPL-2018',
                        real_data = real_moodle_data['2019_demo'])
                        
    hx_list[4] = HX(tube_number = 16,
                        baffle_number = 12,
                        pitch = 10e-3,
                        tube_length = 221.7e-3,
                        plenum_length_1 = 53e-3,
                        plenum_length_2 = 20e-3,
                        baffle_gap = 11.34e-3,
                        baffle_type = 'across_c',
                        tube_layout='t',
                        shell_passes = 2,
                        tube_bundle_diameter= 50e-3,
                        tube_passes = 4,
                        baffle_spacing_in = 40.5e-3,
                        baffle_spacing_out = 15.1e-3,
                        design_year = 2019,
                        pump_year = 2019,
                        T_inh = 57.1,
                        T_inc = 19.6,
                        name = '2019_A1',
                        real_data = real_moodle_data['2019_A1'])
                        #2019 Group A Run 1

    hx_list[5] = HX(tube_number = 14,
                       baffle_number = 11,
                       pitch = 10e-3,
                       tube_length = 226e-3,
                       plenum_length_1 = 55e-3,
                       plenum_length_2 = 20e-3,
                       baffle_gap = 14e-3,
                       baffle_type = 'across_c',
                       tube_layout='t',
                       shell_passes = 2,
                       tube_bundle_diameter= 40e-3,
                       tube_passes = 2,
                       baffle_spacing_in = 40.25e-3,
                       baffle_spacing_out = 15.25e-3,
                       design_year = 2019,
                       pump_year = 2019,
                       T_inh = 54.3,
                       T_inc = 19.3,
                       name = '2019_B1',
                       real_data = real_moodle_data['2019_B1'])
                       #2019 Group B Run 1

    hx_list[6] = HX(tube_number = 16,
                        baffle_number = 6,
                        pitch = 10e-3,
                        tube_length = 222e-3,
                        plenum_length_1 = 53e-3,
                        plenum_length_2 = 24e-3,
                        baffle_gap = 17.34e-3,
                        baffle_type = 'across_c',
                        tube_layout='t',
                        shell_passes = 1,
                        tube_bundle_diameter= 50e-3,
                        tube_passes = 2,
                        baffle_spacing_in = 38.25e-3,
                        baffle_spacing_out = 38.25e-3,
                        design_year = 2019,
                        pump_year = 2019,
                        T_inh = 52.5,
                        T_inc = 19.4,
                        name = '2019_C1',
                        real_data = real_moodle_data['2019_C1'])
                        #2019 Group C Run 1

    hx_list[7] = HX(tube_number = 20,
                        baffle_number = 6,
                        pitch = 11e-3,
                        tube_length = 193e-3,
                        plenum_length_1 = 53e-3,
                        plenum_length_2 = 23e-3,
                        baffle_gap = 12.8e-3,
                        baffle_type = 'across_c',
                        tube_layout='t',
                        shell_passes = 2,
                        tube_bundle_diameter= 54e-3,
                        tube_passes = 2,
                        baffle_spacing_in = 40.75e-3,
                        baffle_spacing_out = 23.25e-3,
                        design_year = 2019,
                        pump_year = 2019,
                        T_inh = 52.7,
                        T_inc = 19.0,
                        name = '2019_D1',
                        real_data = real_moodle_data['2019_D1'])
                        #2019 Group D Run 1


    hx_list[8] = HX(tube_number = 19,
                        baffle_number = 8,
                        pitch = 11e-3,
                        tube_length = 212e-3,
                        plenum_length_1 = 55e-3,
                        plenum_length_2 = 25e-3,
                        baffle_gap = 11.22e-3,
                        baffle_type = 'across_c',
                        tube_layout='t',
                        shell_passes = 1,
                        tube_bundle_diameter= 50e-3,
                        tube_passes = 2,
                        baffle_spacing_in = 35.5e-3,
                        baffle_spacing_out = 35.5e-3,
                        design_year = 2019,
                        pump_year = 2019,
                        T_inh = 53.9,
                        T_inc = 19.4,
                        name = '2019_A2',
                        real_data = real_moodle_data['2019_A2'])
                        #2019 Group A Run 2

    hx_list[9] = HX(tube_number = 14,
                       baffle_number = 11,
                       pitch = 10e-3,
                       tube_length = 226e-3,
                       plenum_length_1 = 55e-3,
                       plenum_length_2 = 20e-3,
                       baffle_gap = 14e-3,
                       baffle_type = 'across_c',
                       tube_layout='t',
                       shell_passes = 2,
                       tube_bundle_diameter= 40e-3,
                       tube_passes = 2,
                       baffle_spacing_in = 40.25e-3,
                       baffle_spacing_out = 15.25e-3,
                       design_year = 2019,
                       pump_year = 2019,
                       T_inh = 55.2,
                       T_inc = 19.5,
                       name = '2019_B2',
                       real_data = real_moodle_data['2019_B2'])
                       #2019 Group B Run 2

    hx_list[10] = HX(tube_number = 16,
                        baffle_number = 6,
                        pitch = 10e-3,
                        tube_length = 222e-3,
                        plenum_length_1 = 53e-3,
                        plenum_length_2 = 24e-3,
                        baffle_gap = 17.34e-3,
                        baffle_type = 'across_c',
                        tube_layout='t',
                        shell_passes = 1,
                        tube_bundle_diameter= 50e-3,
                        tube_passes = 2,
                        baffle_spacing_in = 38.25e-3,
                        baffle_spacing_out = 38.25e-3,
                        design_year = 2019,
                        pump_year = 2019,
                        T_inh = 53.0,
                        T_inc = 19.6,
                        name = '2019_C2',
                        real_data = real_moodle_data['2019_C2'])
                        #2019 Group C Run 2

    hx_list[11] = HX(tube_number = 20,
                        baffle_number = 6,
                        pitch = 11e-3,
                        tube_length = 193e-3,
                        plenum_length_1 = 53e-3,
                        plenum_length_2 = 23e-3,
                        baffle_gap = 6.4e-3,
                        baffle_type = 'across_c',
                        tube_layout='t',
                        shell_passes = 2,
                        tube_bundle_diameter= 54e-3,
                        tube_passes = 2,
                        baffle_spacing_in = 40.75e-3,
                        baffle_spacing_out = 23.25e-3,
                        design_year = 2019,
                        pump_year = 2019,
                        T_inh = 53.2,
                        T_inc = 19.4,
                        name = '2019_D2',
                        real_data = real_moodle_data['2019_D2'])
                        #2019 Group D Run 2

    hx_list[12] = HX(tube_number=20,
                    baffle_number=8,
                    pitch=13e-3,
                    tube_length=141e-3,
                    plenum_length_1=53e-3,
                    plenum_length_2=24e-3,
                    baffle_gap=18.92e-3,
                    baffle_type='across_c',
                    tube_layout='t',
                    shell_passes=2,
                    tube_bundle_diameter=58e-3,
                    tube_passes=4,
                    design_year=2018,
                    pump_year=2018,
                    T_inh=53.5,
                    T_inc=19.6,
                    baffle_spacing_in=26.7e-3,
                    baffle_spacing_out=26.7e-3,
                    name='2018_A1',
                    real_data=real_moodle_data['2018_A1']
                    )

    hx_list[13] = HX(tube_number=24,
                    baffle_number=6,
                    pitch=10e-3,
                    tube_length=130e-3,
                    plenum_length_1=53e-3,
                    plenum_length_2=24e-3,
                    baffle_gap=18.75e-3,
                    baffle_type='across_c',
                    tube_layout='t',
                    shell_passes=2,
                    tube_bundle_diameter=43.5e-3,
                    tube_passes=2,
                    design_year=2018,
                    pump_year=2018,
                    T_inh=48.4,
                    T_inc=22.4,
                    baffle_spacing_in=45e-3,
                    baffle_spacing_out=45e-3,
                    name='2018_B1',
                    real_data=real_moodle_data['2018_B1']
                    )

    hx_list[14] = HX(tube_number=20,
                    baffle_number=5,
                    pitch=10e-3,
                    tube_length=156e-3,
                    plenum_length_1=53e-3,
                    plenum_length_2=24e-3,
                    baffle_gap=15e-3,
                    baffle_type='across_c',
                    tube_layout='t',
                    shell_passes=2,
                    tube_bundle_diameter=57e-3,
                    tube_passes=2,
                    design_year=2018,
                    pump_year=2018,
                    T_inh=52.5,
                    T_inc=20.5,
                    baffle_spacing_in=42e-3,
                    baffle_spacing_out=42e-3,
                    name='2018_C1',
                    real_data=real_moodle_data['2018_C1']
                    )

    hx_list[15] = HX(tube_number=20,
                    baffle_number=5,
                    pitch=10e-3,
                    tube_length=156e-3,
                    plenum_length_1=53e-3,
                    plenum_length_2=24e-3,
                    baffle_gap=15e-3,
                    baffle_type='across_c',
                    tube_layout='t',
                    shell_passes=2,
                    tube_bundle_diameter=57e-3,
                    tube_passes=2,
                    design_year=2018,
                    pump_year=2018,
                    T_inh=51.3,
                    T_inc=21.3,
                    baffle_spacing_in=42e-3,
                    baffle_spacing_out=42e-3,
                    name='2018_C2',
                    real_data=real_moodle_data['2018_C2']
                    )

    hx_list[16] = HX(tube_number=20,
                    baffle_number=8,
                    pitch=13e-3,
                    tube_length=141e-3,
                    plenum_length_1=53e-3,
                    plenum_length_2=24e-3,
                    baffle_gap=18.92e-3,
                    baffle_type='across_c',
                    tube_layout='t',
                    shell_passes=2,
                    tube_bundle_diameter=58e-3,
                    tube_passes=4,
                    design_year=2018,
                    pump_year=2018,
                    T_inh=50.0,
                    T_inc=21.8,
                    baffle_spacing_in=26.7e-3,
                    baffle_spacing_out=26.7e-3,
                    name='2018_A2',
                    real_data=real_moodle_data['2018_A2']
                    )

    # hx_list[8] = HX(tube_number = 19,
    #                     baffle_number = 8,
    #                     pitch = 11e-3,
    #                     tube_length = 212e-3,
    #                     plenum_length_1 = 55e-3,
    #                     plenum_length_2 = 25e-3,
    #                     baffle_gap = 11.22e-3,
    #                     baffle_type = 'across_c',
    #                     tube_layout='t',
    #                     shell_passes = 1,
    #                     tube_bundle_diameter= 50e-3,
    #                     tube_passes = 2,
    #                     baffle_spacing_in = 35.5e-3,
    #                     baffle_spacing_out = 35.5e-3,
    #                     design_year = 2019,
    #                     pump_year = 2019,
    #                     T_inh = 53.4,
    #                     T_inc = 18.9,
    #                     name = '2019_E1',
    #                     real_data = real_moodle_data['2019_E1'])
    #                     #2019 Group E Run 1

    # hx_list[13] = HX(tube_number = 19,
    #                     baffle_number = 8,
    #                     pitch = 11e-3,
    #                     tube_length = 212e-3,
    #                     plenum_length_1 = 55e-3,
    #                     plenum_length_2 = 25e-3,
    #                     baffle_gap = 11.22e-3,
    #                     baffle_type = 'across_c',
    #                     tube_layout='t',
    #                     shell_passes = 1,
    #                     tube_bundle_diameter= 50e-3,
    #                     tube_passes = 2,
    #                     baffle_spacing_in = 35.5e-3,
    #                     baffle_spacing_out = 35.5e-3,
    #                     design_year = 2019,
    #                     pump_year = 2019,
    #                     T_inh = 53.8,
    #                     T_inc = 19.1,
    #                     name = '2019_E2',
    #                     real_data = real_moodle_data['2019_E2'])
    #                     #2019 Group E Run 2

    # hx_list[14] = HX(tube_number=8,
    #             baffle_number=12,
    #             pitch=10e-3,
    #             tube_length=261e-3,
    #             plenum_length_1=41e-3,
    #             plenum_length_2=41e-3,
    #             baffle_gap=15.75e-3,
    #             baffle_type='across_c',
    #             tube_layout='t',
    #             shell_passes=1,
    #             tube_bundle_diameter=45e-3,
    #             tube_passes=1,
    #             design_year=2017,
    #             pump_year=2017,
    #             T_inh=60,
    #             T_inc=20,
    #             baffle_spacing_in=34e-3,
    #             baffle_spacing_out=34e-3,
    #             name='2017_A',
    #             real_data=None
    #             )

    # hx_list[15] = HX(tube_number=14,
    #                 baffle_number=12,
    #                 pitch=15e-3,
    #                 tube_length=235e-3.3,
    #                 plenum_length_1=41e-3,
    #                 plenum_length_2=23e-3,
    #                 baffle_gap=25.6e-3,
    #                 baffle_type='across_c',
    #                 tube_layout='t',
    #                 shell_passes=2,
    #                 tube_bundle_diameter=55e-3,
    #                 tube_passes=2,
    #                 design_year=2017,
    #                 pump_year=2017,
    #                 T_inh=60,
    #                 T_inc=20,
    #                 baffle_spacing_in=33e-3,
    #                 baffle_spacing_out=16e-3,
    #                 name='2017_B',
    #                 real_data=None
    #                 )

    # hx_list[16] = HX(tube_number=14,
    #                 baffle_number=8,
    #                 pitch=10e-3,
    #                 tube_length=192e-3,
    #                 plenum_length_1=41e-3,
    #                 plenum_length_2=24e-3,
    #                 baffle_gap=15.75e-3,
    #                 baffle_type='across_c',
    #                 tube_layout='t',
    #                 shell_passes=1,
    #                 tube_bundle_diameter=42e-3,
    #                 tube_passes=2,
    #                 design_year=2017,
    #                 pump_year=2017,
    #                 T_inh=60,
    #                 T_inc=20,
    #                 baffle_spacing_in=32e-3,
    #                 baffle_spacing_out=32e-3,
    #                 name='2017_C',
    #                 real_data=None
    #                 )

    # hx_list[17] = HX(tube_number=20,
    #                 baffle_number=6,
    #                 pitch=10e-3,
    #                 tube_length=147e-3,
    #                 plenum_length_1=41e-3,
    #                 plenum_length_2=24e-3,
    #                 baffle_gap=15.75e-3,
    #                 baffle_type='across_c',
    #                 tube_layout='t',
    #                 shell_passes=1,
    #                 tube_bundle_diameter=58e-3,
    #                 tube_passes=2,
    #                 design_year=2017,
    #                 pump_year=2017,
    #                 T_inh=60,
    #                 T_inc=20,
    #                 baffle_spacing_in=32e-3,
    #                 baffle_spacing_out=32e-3,
    #                 name='2017_D',
    #                 real_data=None
    #                 )

    # hx_list[18] = HX(tube_number=24,
    #                 baffle_number=7,
    #                 pitch=10e-3,
    #                 tube_length=132e-3,
    #                 plenum_length_1=31e-3,
    #                 plenum_length_2=24e-3,
    #                 baffle_gap=15.8e-3,
    #                 baffle_type='across_c',
    #                 tube_layout='t',
    #                 shell_passes=1,
    #                 tube_bundle_diameter=60e-3,
    #                 tube_passes=2,
    #                 design_year=2017,
    #                 pump_year=2017,
    #                 T_inh=60,
    #                 T_inc=20,
    #                 baffle_spacing_in=32e-3,
    #                 baffle_spacing_out=32e-3,
    #                 name='2017_E',
    #                 real_data=None
    #                 )




    if heat_exchanger == None:
        return hx_list
    else:
        hx_singular = {}
        hx_singular[heat_exchanger] = hx_list[heat_exchanger]
        return hx_singular



def moodle_performance_dict(heat_exchanger = None):

    p_data = {}

    p_data['hx_y2018_p2022'] = {'Name':'real_data',
                                'T cold out (C)':26.1,   
                                'T hot out (C)':46.8,
                                'mdot_cold (l/s)':0.454,
                                'mdot_hot (l/s)':0.472,
                                'dP_cold (bar)':0.277,
                                'dP_hot (bar)':0.137,
                                'Q_NTU (kW)':12.460,
                                'Q_NTU_corr (kW)':14.79,
                                'eff_NTU':0.197,
                                'mass (kg)':1.466
                                }

    p_data['hx_y2017B_p2022'] = {'Name':'real_data',
                                'T cold out (C)':30.1,   
                                'T hot out (C)':43.4,
                                'mdot_cold (l/s)':0.250,
                                'mdot_hot (l/s)':0.425,
                                'dP_cold (bar)':0.471,
                                'dP_hot (bar)':0.212,
                                'Q_NTU (kW)':8.520,
                                'Q_NTU_corr (kW)':13.01,
                                'eff_NTU':0.314,
                                'mass (kg)':1.021
                                }

    p_data['2019_demo'] = {'Name':'real_data',
                        'T cold out (C)':21,
                       'T hot out (C)':45.3,
                       'mdot_cold (l/s)':0.525,
                       'mdot_hot (l/s)':0.490,
                        'dP_cold (bar)':0.356,
                       'dP_hot (bar)':0.163,
                       'Q_NTU (kW)':15.82,
                       'Q_NTU_corr (kW)':16.06,
                       'eff_NTU':0.2,
                       'mass (kg)':1.466
                       }

    p_data['2019_A1'] = {'Name':'real_data',
                        'T cold out (C)':26.4,
                        'T hot out (C)':47.2,
                        'mdot_cold (l/s)':0.417,
                        'mdot_hot (l/s)':0.299,
                        'dP_cold (bar)':0.437,
                        'dP_hot (bar)':0.391,
                        'Q_NTU (kW)':11.98,
                        'Q_NTU_corr (kW)':12.78,
                        'eff_NTU':0.26,
                        'mass (kg)':1.15
                        }

    p_data['2019_B1'] = {'Name':'real_data',
                        'T cold out (C)':26,
                        'T hot out (C)':48,
                        'mdot_cold (l/s)':0.375,
                        'mdot_hot (l/s)':0.462,
                        'dP_cold (bar)':0.458,
                        'dP_hot (bar)':0.219,
                        'Q_NTU (kW)':11.22,
                        'Q_NTU_corr (kW)':12.82,
                        'eff_NTU':0.21,
                        'mass (kg)':1.03
                        }

    p_data['2019_C1'] = {'Name':'real_data',
                        'T cold out (C)':24.1,
                        'T hot out (C)':46.2,
                        'mdot_cold (l/s)':0.592,
                        'mdot_hot (l/s)':0.462,
                        'dP_cold (bar)':0.275,
                        'dP_hot (bar)':0.212,
                        'Q_NTU (kW)':11.77,
                        'Q_NTU_corr (kW)':14.23,
                        'eff_NTU':0.19,
                        'mass (kg)':1.09
                        }

    p_data['2019_D1'] = {'Name':'real_data',
                        'T cold out (C)':25,
                        'T hot out (C)':43.6,
                        'mdot_cold (l/s)':0.575,
                        'mdot_hot (l/s)':0.349,
                        'dP_cold (bar)':0.362,
                        'dP_hot (bar)':0.348,
                        'Q_NTU (kW)':13.71,
                        'Q_NTU_corr (kW)':16.27,
                        'eff_NTU':0.28,
                        'mass (kg)':1.16
                        }

    p_data['2019_E1'] = {'Name':'real_data',
                        'T cold out (C)':25,
                        'T hot out (C)':46.2,
                        'mdot_cold (l/s)':0.608,
                        'mdot_hot (l/s)':0.483,
                        'dP_cold (bar)':0.325,
                        'dP_hot (bar)':0.2,
                        'Q_NTU (kW)':14.87,
                        'Q_NTU_corr (kW)':17.24,
                        'eff_NTU':0.22,
                        'mass (kg)':1.11
                        }

    p_data['2019_E2'] = {'Name':'real_data',
                        'T cold out (C)':25.4,
                        'T hot out (C)':46.6,
                        'mdot_cold (l/s)':0.608,
                        'mdot_hot (l/s)':0.479,
                        'dP_cold (bar)':0.322,
                        'dP_hot (bar)':0.201,
                        'Q_NTU (kW)':15.07,
                        'Q_NTU_corr (kW)':17.37,
                        'eff_NTU':0.22,
                        'mass (kg)':1.11
                        }

    p_data['2019_D2'] = {'Name':'real_data',
                        'T cold out (C)':25.2,
                        'T hot out (C)':44.6,
                        'mdot_cold (l/s)':0.575,
                        'mdot_hot (l/s)':0.354,
                        'dP_cold (bar)':0.359,
                        'dP_hot (bar)':0.348,
                        'Q_NTU (kW)':13.2,
                        'Q_NTU_corr (kW)':15.62,
                        'eff_NTU':0.27,
                        'mass (kg)':1.16
                        }

    p_data['2019_C2'] = {'Name':'real_data',
                        'T cold out (C)':24.2,
                        'T hot out (C)':46.2,
                        'mdot_cold (l/s)':0.592,
                        'mdot_hot (l/s)':0.476,
                        'dP_cold (bar)':0.269,
                        'dP_hot (bar)':0.199,
                        'Q_NTU (kW)':12.32,
                        'Q_NTU_corr (kW)':14.76,
                        'eff_NTU':0.19,
                        'mass (kg)':1.09
                        }

    p_data['2019_B2'] = {'Name':'real_data',
                        'T cold out (C)':26.1,
                        'T hot out (C)':48.8,
                        'mdot_cold (l/s)':0.375,
                        'mdot_hot (l/s)':0.462,
                        'dP_cold (bar)':0.449,
                        'dP_hot (bar)':0.223,
                        'Q_NTU (kW)':11.23,
                        'Q_NTU_corr (kW)':12.59,
                        'eff_NTU':0.2,
                        'mass (kg)':1.03
                        }

    p_data['2019_A2'] = {'Name':'real_data',
                        'T cold out (C)':25.4,
                        'T hot out (C)':44.9,
                        'mdot_cold (l/s)':0.417,
                        'mdot_hot (l/s)':0.299,
                        'dP_cold (bar)':0.422,
                        'dP_hot (bar)':0.4,
                        'Q_NTU (kW)':10.73,
                        'Q_NTU_corr (kW)':12.44,
                        'eff_NTU':0.25,
                        'mass (kg)':1.15
                        }

    p_data['2018_A1'] = {'Name':'real_data',
                        'T cold out (C)':27.2,
                        'T hot out (C)':48.3,
                        'mdot_cold (l/s)':0.317,
                        'mdot_hot (l/s)':0.365,
                        'dP_cold (bar)':0.202,
                        'dP_hot (bar)':0.292,
                        'Q_NTU (kW)':8.91,
                        'Q_NTU_corr (kW)':10.52,
                        'eff_NTU':0.2,
                        'mass (kg)':0.898
                        }

    p_data['2018_C1'] = {'Name':'real_data',
                        'T cold out (C)':28.5,
                        'T hot out (C)':47.8,
                        'mdot_cold (l/s)':0.32,
                        'mdot_hot (l/s)':0.421,
                        'dP_cold (bar)':0.241,
                        'dP_hot (bar)':0.179,
                        'Q_NTU (kW)':9.39,
                        'Q_NTU_corr (kW)':11.74,
                        'eff_NTU':0.22,
                        'mass (kg)':0.894
                        }

    p_data['2018_C2'] = {'Name':'real_data',
                        'T cold out (C)':28.6,
                        'T hot out (C)':47,
                        'mdot_cold (l/s)':0.323,
                        'mdot_hot (l/s)':0.421,
                        'dP_cold (bar)':0.244,
                        'dP_hot (bar)':0.179,
                        'Q_NTU (kW)':8.62,
                        'Q_NTU_corr (kW)':11.5,
                        'eff_NTU':0.22,
                        'mass (kg)':0.894
                        }

    p_data['2018_A2'] = {'Name':'real_data',
                        'T cold out (C)':27.9,
                        'T hot out (C)':44.7,
                        'mdot_cold (l/s)':0.349,
                        'mdot_hot (l/s)':0.349,
                        'dP_cold (bar)':0.199,
                        'dP_hot (bar)':0.289,
                        'Q_NTU (kW)':8.22,
                        'Q_NTU_corr (kW)':11.66,
                        'eff_NTU':0.2,
                        'mass (kg)':0.898
                        }

    p_data['2018_B1'] = {'Name':'real_data',
                        'T cold out (C)':30.3,
                        'T hot out (C)':45.5,
                        'mdot_cold (l/s)':0.292,
                        'mdot_hot (l/s)':0.464,
                        'dP_cold (bar)':0.285,
                        'dP_hot (bar)':0.13,
                        'Q_NTU (kW)':7.55,
                        'Q_NTU_corr (kW)':11.61,
                        'eff_NTU':0.24,
                        'mass (kg)':0.935
                        }


    # #correct


    # #correct


 

    if heat_exchanger == None:
        return p_data
    else:
        hx_singular = {}
        hx_singular[heat_exchanger] = p_data[heat_exchanger]
        return hx_singular



def dict_2022(heat_exchanger = None):

    hx_list = {}



    hx_list['9'] = HX(tube_number=14,
                    baffle_number=12,
                    pitch=14e-3,
                    tube_length=236e-3,
                    plenum_length_1=51e-3,
                    plenum_length_2=51e-3,
                    baffle_gap=12.00e-3,
                    baffle_type='across_c',
                    tube_layout='t',
                    shell_passes=1,
                    tube_bundle_diameter=56e-3,
                    tube_passes=1,
                    design_year=2022,
                    pump_year=2022,
                    T_inh=60,
                    T_inc=20,
                    baffle_spacing_in=39.15e-3,
                    baffle_spacing_out=39.15e-3,
                    name='2022_A',
                    real_data=None
                    )

    hx_list['10'] = HX(tube_number=16,
                     baffle_number=6,
                     pitch=10e-3,
                     tube_length=219e-3,
                     plenum_length_1=50e-3,
                     plenum_length_2=50e-3,
                     baffle_gap=17.22e-3,
                     baffle_type='across_c',
                     tube_layout='t',
                     shell_passes=1,
                     tube_bundle_diameter=52e-3,
                     tube_passes=1,
                     design_year=2022,
                     pump_year=2022,
                     T_inh=60,
                     T_inc=20,
                     baffle_spacing_in=39.5e-3,
                     baffle_spacing_out=39.5e-3,
                     name='2022_B',
                     real_data=None
                     )

    hx_list['11'] = HX(tube_number=20,
                     baffle_number=4,
                     pitch=10e-3,
                     tube_length=174e-3,
                     plenum_length_1=50e-3,
                     plenum_length_2=46.5e-3,
                     baffle_gap=9.75e-3,
                     baffle_type='across_c',
                     tube_layout='t',
                     shell_passes=1,
                     tube_bundle_diameter=60e-3,
                     tube_passes=2,
                     design_year=2022,
                     pump_year=2022,
                     T_inh=60,
                     T_inc=20,
                     baffle_spacing_in=39e-3,
                     baffle_spacing_out=39e-3,
                     name='2022_C',
                     real_data=None)

    hx_list['12'] = HX(tube_number=12,
                     baffle_number=8,
                     pitch=12e-3,
                     tube_length=259e-3,
                     plenum_length_1=53e-3,
                     plenum_length_2=30e-3,
                     baffle_gap=15.36e-3,
                     baffle_type='across_c',
                     tube_layout='t',
                     shell_passes=1,
                     tube_bundle_diameter=52e-3,
                     tube_passes=2,
                     design_year=2022,
                     pump_year=2022,
                     T_inh=60,
                     T_inc=20,
                     baffle_spacing_in=49e-3,
                     baffle_spacing_out=49e-3,
                     name='2022_D',
                     real_data=None
                     )

    hx_list['13'] = HX(tube_number=24,
                     baffle_number=4,
                     pitch=10e-3,
                     tube_length=145e-3,
                     plenum_length_1=53e-3,
                     plenum_length_2=41e-3,
                     baffle_gap=16.30e-3,
                     baffle_type='across_c',
                     tube_layout='t',
                     shell_passes=1,
                     tube_bundle_diameter=62.5e-3,
                     tube_passes=4,
                     design_year=2022,
                     pump_year=2022,
                     T_inh=60,
                     T_inc=20,
                     baffle_spacing_in=40e-3,
                     baffle_spacing_out=40e-3,
                     name='2022_E',
                     real_data=None
                     )

    hx_list['14'] = HX(tube_number=20,
                     baffle_number=6,
                     pitch=10e-3,
                     tube_length=172e-3,
                     plenum_length_1=58e-3,
                     plenum_length_2=23e-3,
                     baffle_gap=20.00e-3,
                     baffle_type='across_c',
                     tube_layout='t',
                     shell_passes=2,
                     tube_bundle_diameter=55e-3,
                     tube_passes=4,
                     design_year=2022,
                     pump_year=2022,
                     T_inh=60,
                     T_inc=20,
                     baffle_spacing_in=41.5e-3,
                     baffle_spacing_out=21.5e-3,
                     name='2022_F',
                     real_data=None
                     )
    
    hx_list['15'] = HX(tube_number = 13,
                        baffle_number = 8,
                        pitch = 12e-3,
                        tube_length = 212e-3,
                        plenum_length_1 = 53.5e-3,
                        plenum_length_2 = 53.5e-3,
                        baffle_gap = 13.9e-3,
                        baffle_type = 'across_c',
                        tube_layout='t',
                        shell_passes = 1,
                        tube_bundle_diameter= 56e-3,
                        tube_passes = 1,
                        baffle_spacing_in = 38e-3,
                        baffle_spacing_out = 38e-3,
                        design_year = 2019,
                        pump_year = 2019,
                        T_inh = 53,
                        T_inc = 19.5,
                        name = 'BLIND_TEST',
                        real_data = None
                        )
                        #Blind test case


    if heat_exchanger == None:
        return hx_list
    else:
        hx_singular = {}
        hx_singular[heat_exchanger] = hx_list[heat_exchanger]
        return hx_singular

#endregion

#OTHER TOOLS
#region

# HX_MOODLE_DATA: tool which outputs a dataframe of data for the moodle heat exchangers

def predict_hx(data = 'moodle',heat_exchanger=None,k_array = np.array([1,1,1,1,1,1,1,1,1,1]),sort_data = False):

    if data == 'moodle':
        hx_list = heat_exchanger_dict(heat_exchanger)
    elif data == '2022':
        hx_list = dict_2022(heat_exchanger)
    else:
        print('choose a correct dataset')
    hx_data = pd.DataFrame()

    for hxi in hx_list:
        hx = hx_list[hxi]

        if data == 'moodle':

            real_data = hx.real_data
            hx_data = hx_data.append(real_data, ignore_index = True)

        performance = hx_design(hx,k_array) 
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
                'Q_NTU_corr (kW)',
                'eff_NTU',
                'mass (kg)'
                    ]]
    if sort_data == True:
        hx_data = hx_data.sort_values(by="Q_NTU_corr (kW)", ascending=False).head()
    with pd.option_context('display.max_rows', None, 'display.max_columns', None,"display.precision", 3):  # more options can be specified also
        print(hx_data)

    hx_data.to_excel(f"prediction_data_{data}final5.xlsx", sheet_name="prediction_data", index=False)



def optimiser(n = 10):
    #brute force optimisation with a few checks to eliminate cases as early as possible

    baffle_type = 'across_c' #assume this is best
    tube_layout = 't'        #assume this is best

    hx_designs = {}
    hx_data = pd.DataFrame()
    
    design_no = 0
    run1 = 'run_l'

    tp_array = np.array([4])
    sp_array = np.array([2])
    pl1_array = np.array([55e-3]) #np.linspace(41e-3,100e-3,n)
    bsi_array = np.array([40e-3]) #np.linspace(41e-3,50e-3,n)
    bn_array = np.array([6,7,8,9,10,11,12,13,14])   #PLEASE DONT HAVE ZERO OR 1 BAFFLE  
    bg_array = np.array([10e-3,15e-3,20e-3]) #15e-3,20e-3,
    p_array = np.array([10e-3])
    pl2_array = np.array([20e-3])
    bso_array = np.array([20e-3])

    packing_density = np.pi/np.sqrt(12)
    r = 64e-3/2
    segment_area = r**2 * np.arccos(0.8660) - (0.8660*r)*(2*r*0.13397*r - (0.13397*r)**2)**0.5
    max_tube_area = (np.pi*(r**2) - 6*segment_area)*(packing_density - 0.1)

    initial_length = 85e-3
    tube = pipe(6e-3,8e-3,0.20,3.5,initial_length)

    for pitch in p_array:
        print('pitch: ',pitch)

        eff_tube_area = np.pi*(pitch/2)**2
        max_no_tubes_from_area = int(max_tube_area/eff_tube_area)
        max_no_tube_from_mass = int(1.1/tube.mass)
        max_no_tubes = min(max_no_tubes_from_area,max_no_tube_from_mass)
        max_no_tubes = 20
        min_no_tubes = int(max_no_tubes*0.4)
        tn_array = np.array(range(int(min_no_tubes),int(max_no_tubes + 1)))
        tn_array = np.array([20])

        for tube_passes in tp_array:
            #print('tube_passes: ',tube_passes)

            # if tube_passes%2 == 0:
            #     l_min = 10e-3
            # else:
            #     l_min = 41e-3
             #np.linspace(l_min,100e-3,n)  20e-3,30e-3,
        
            for shell_passes in sp_array:
                #print('shell_passes: ',shell_passes)

                # if shell_passes%2 == 0:
                #     bso_min = 10e-3
                # else:
                #     bso_min = 41e-3
                 #np.linspace(bso_min,50e-3,n)  20e-3,30e-3,
                for tube_number in tn_array:
                    #print('tube_number',tube_number)

                    for plenum_length_1 in pl1_array:                        

                        for plenum_length_2 in pl2_array:

                            for baffle_spacing_in in bsi_array:

                                for baffle_spacing_out in bso_array:
                    
                                    max_tube_length = (0.35 - plenum_length_1 - plenum_length_2)
                                    min_tube_length = baffle_spacing_in + baffle_spacing_out + 0.01
                                    tl_array = np.linspace(min_tube_length,max_tube_length,n)
                                    n_total = tp_array.size*sp_array.size*tl_array.size*pl1_array.size*bsi_array.size*bn_array.size*bg_array.size*p_array.size*bso_array.size*tn_array.size*pl2_array.size

                                    for tube_length in tl_array:

                                        tube.l = tube_length

                                        if (tube_number * (tube_length+0.003) <= 3.5) and (tube_number * tube.mass <= 1.1):               
                                        

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
                                                                        name = f'design number = {design_no}'
                                                                        )
                                                    
                                                    if heat_exchanger.total_mass() <= 1.1:
                                                        hx_designs[f'design {design_no}'] = heat_exchanger
                                                        performance = hx_design(heat_exchanger) #,invalid_hx_flag
                                                        design = vars(heat_exchanger)
                                                        performance.update(design)
                                                        #if invalid_hx_flag == False:
                                                        hx_data = hx_data.append(performance, design) 
                                                    

    #order columns nicely
    hx_data = hx_data.sort_values(by="Q_NTU (kW)", ascending=False)[['Name',
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

#endregion

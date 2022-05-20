from re import M
from tkinter import N
from scipy.optimize import fsolve
import numpy as np
from scipy.interpolate import interp1d


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

def Re(V,d,liquid):
    Re = (V*d*liquid.rho)/liquid.mu
    return Re

def hi(V_tube,di,liquid):

    R = Re(V_tube,di,liquid)

    Nu = 0.023*R**0.8*liquid.Pr**0.3

    hi = Nu*liquid.k/di

    return hi

def ho(Cc, m_c, pitch, d_outer, liquid, baffle_spacing, d_shell, d_otl): # bell delaware method

    #from correlations at https://reader.elsevier.com/reader/sd/pii/B9780123735881500097?token=5F15AE11E5C1F5036E719B79796AE27989BCE35E4CC9A873B8C8DF9B21C28E1F52DA6116E24B9C5E327EBC98518B8F9D&originRegion=eu-west-1&originCreation=20220519125049
    a1 = 0.321 
    a2 = -0.388
    a3 = 1.450
    a4 = 0.519
    #for 30 degree layout, i.e. triangular

    #areas
    Sm = baffle_spacing * ((d_shell - d_otl) + ((d_otl - d_outer)*(pitch - d_outer))/pitch) #incorporate number of passes???
    Sb = baffle_spacing * (d_shell - d_otl)
    #Ssb = 0
    #Stb = 0

    G = m_c / Sm
    a = a3 / (1 + 0.14 * ((d_outer * G)/liquid.mu)**a4)
    j = a1 * (1.33 /(pitch/d_outer))**a * ((d_outer * G)/liquid.mu)**a2 #dimensionless quanitity
    phi = 1 #viscosity correction factor

    ho_ideal = j * liquid.cp * G * phi / liquid.Pr**(2/3)
    
    #rs = Ssb / (Ssb + Stb)
    #rl = (Ssb + Stb) / Sm

    Jc = 1#0.55 + 0.72 * Fc #correction factor for baffle window flow
    Jl = 1#0.44 * (1- rs) + (1 - 0.44 * (1 -rs))*np.exp(-2.2 * rl) #correction factor for baffle leakage effects
    Jb = np.exp(-1.35 * (Sb/Sm)) #correction factor for bundle bypass effects
    Jr = 1 #laminar flow correction factor
    Js = 1 #correction factor for unequal baffle spacing

    ho = ho_ideal * Jc * Jl * Jb * Jr * Js
    return ho

def ho1(do, liquid, pitch, baffle_area, A_shell, d_shell, tube_length, baffle_spacing, m_dot):
    #relation from https://reader.elsevier.com/reader/sd/pii/0017931063900371?token=67FD7B1EA2E8B9D710BA94CFEC6DA5F06A07655CAAFAFD0A8CC482CCE90F3D4CC6E9878A64AE9F808EAC4F0844192201&originRegion=eu-west-1&originCreation=20220517160904
    L1 = tube_length  #distance between end plates, not quite 
    L3 = baffle_spacing #distance between baffles
    L2 = L1 - 2*L3 #'length between end baffles' not quite sure how to interpret this

    Sw = A_shell - baffle_area #Sw is the free area for flow to go through in the plane of the baffles
    Se = 0 #Se is the leakage area around baffles, assume zero for now

    Sm = d_shell - 5 * do #maximum free area for flow in cross flow zone MODIFY 5
    Sp = d_shell - 5 * do #minimum free area for flow in cross flow zone

    P = ((pitch - do)/pitch * (do/d_shell))

    Gav = 1/3 * (m_dot/Sw + m_dot/Sm + m_dot/Sp)
    #print(Gav)

    R = Re(Gav,do,liquid)
    #print(R)

    S = (Sw / (Sw + Se))

    F = (L2 + (L1 - L2)*(2*L3/(L1-L2))**0.6)/L1


    Nu = 1.9 * R**0.6 * liquid.Pr**0.3 * P**0.4 * S**2 * F

    ho = Nu*liquid.k/do

    return ho

def ho2(V_shell,do,liquid,tube_layout): #handout relation
    if tube_layout == 't':
        c = 0.2
    elif tube_layout == 's':
        c = 0.15
    else:
        c= 0.2
        print('error, invalid tube layout')

    R = Re(V_shell,do,liquid)

    Nu = c*R**0.6*liquid.Pr**0.3

    ho = Nu*liquid.k/do

    return ho

def LMTD(T_1in,T_2in,T_1out,T_2out):
    lmtd = ((T_2in - T_1out) - (T_2out - T_1in)) / np.log((T_2in - T_1out) / (T_2out - T_1in))
    return lmtd

def U_inside(hi,ho,di,do,L,k_copper = 398):

    ro = do/2
    ri = di/2
    Ao = np.pi*do*L
    Ai = np.pi*di*L

    inv_U = 1/hi + (Ai*np.log(ro/ri))/(2*np.pi*k_copper*L) + Ai/(Ao*ho)
    
    return 1/inv_U

def dP_tube(L,di,liquid,V):

    f = friction_factor(Re(V,di,liquid))

    dP = f*(L/di)*0.5*liquid.rho*V**2
    return dP

def dP_inout(liquid,V,sigma):

    Kc,Ke = KcKe(sigma)

    dP = 0.5*liquid.rho*(V**2)*(Kc + Ke)

    return dP

def dP_nozzle(V,liquid):
    dP = 0.5*liquid.rho*(V**2)
    return dP

def dP_shell(liquid,d_outer,N,tube_layout,d_shell, pitch, baffle_spacing, d_otl, m_c, baffle_cut, baffle_number, baffle_area):

    if tube_layout == 't':
        theta = 30
        b1 = 0.372
        b2 = -0.123
        b3 = 7
        b4 = 0.5
    if tube_layout == 's':
        theta = 0
        b1 = 0.0815
        b2 = 0.022
        b3 = 6.30
        b4 = 0.378

    Sm = baffle_spacing * ((d_shell - d_otl) + ((d_otl - d_outer)*(pitch - d_outer))/pitch)
    Swg = (1/4) * np.pi * d_shell **2 - baffle_area 
    Sw = Swg #- N * (Swg/((1/4) * np.pi * d_shell **2)) * d_outer**2 * np.pi * 0.25
    Sb = baffle_spacing * (d_shell - d_otl)

    #ideal pressure drop
    G = m_c/Sm
    b = b3 / (1 + 0.14 * ((d_outer * G)/liquid.mu)**b4)
    Nc = d_shell * (1 - 2 * baffle_cut) / (pitch * np.cos(theta))
    phi = 1
    f = b1 * (1.33/(pitch/d_outer))**b * ((d_outer * G)/liquid.mu)**b2
    dp_ideal = (2 * f * Nc * G**2)/(1 * liquid.rho * phi)
    print(dp_ideal)

    Rs = 1 #uniform baffle spacing
    Rb = np.exp(-3.7*(Sb/Sm))
    Rl = 1 #assume no leakage
    nb = baffle_number
    Ncw = (0.8 * baffle_cut * d_shell) / (pitch * np.cos(theta))

    dpw_ideal = ((2 + 0.6 * Ncw) * m_c**2) / (2 * 1 * liquid.rho * Sm * Sw)
    print(dpw_ideal)

    dP = ((nb -1)*dp_ideal*Rb + nb*dpw_ideal)*Rl + 2*dp_ideal*(1 + Ncw/Nc)*Rb*Rs

    return dP

def dP_shell_0(V,liquid,do,N,tube_layout):

    Re_shell = Re(V,do,liquid)

    if tube_layout == 't':
        a = 0.2
    elif tube_layout == 's':
        a = 0.34
    else:
        a = 0.2
        print('error, invalid tube layout')

    dP = 4*a*Re_shell**(-0.15)*N*liquid.rho*(V**2)

    return dP

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

def mdot_dP(m_dot,dP_ovr,side,liquid,year):

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

    if year == 2019:
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

    if year == 2018:
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

    if year == 2017:
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

    mdot_from_dP = interp1d(mdot_dP_array[:,1],mdot_dP_array[:,0],fill_value='extrapolate',kind = 'cubic')
    dP_from_mdot = interp1d(mdot_dP_array[:,0],mdot_dP_array[:,1],fill_value='extrapolate',kind = 'cubic')

    dP_new = dP_from_mdot(m_dot/(liquid.rho/1000))
    m_dot_new = mdot_from_dP(dP_ovr)*(liquid.rho/1000)

    rel_e_dP = (dP_new - dP_ovr)/dP_ovr
    rel_e_mdot = (m_dot_new - m_dot)/m_dot

    m_dot = (m_dot + m_dot_new)/2

    return m_dot,rel_e_dP,rel_e_mdot

def F(T_inc,T_inh,T_outc,T_outh,shell_passes):

    N = shell_passes
    #N = number of shell passes. 2M = tube passes per shell

    P = (T_outh - T_inh)/(T_inc - T_inh)

    R = (T_inc - T_outc)/(T_outh - T_inh)

    S = ((R**2 + 1)**0.5)/(R - 1)

    W = ((1-P*R)/(1 - P))**(1/N)

    F = (S*np.log(W))/np.log((1 + W - S + S*W)/(1 + W + S - S*W))

    return F

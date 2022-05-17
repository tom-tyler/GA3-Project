from re import M
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

def ho(do, tube_number, liquid, d_shell, pitch, baffle_area, A_shell, tube_length, thickness, baffle_number, m_dot):
    #relation from https://reader.elsevier.com/reader/sd/pii/0017931063900371?token=67FD7B1EA2E8B9D710BA94CFEC6DA5F06A07655CAAFAFD0A8CC482CCE90F3D4CC6E9878A64AE9F808EAC4F0844192201&originRegion=eu-west-1&originCreation=20220517160904
    L1 = tube_length - 2*thickness #distance between end plates
    L3 = L1/baffle_number - thickness #distance between baffles
    L2 = L1 - 2*L3 - 2*thickness#'length between end baffles' not quite sure how to interpret this

    Sw = A_shell - baffle_area #Sw is the free area for flow to go through in the plane of the baffles
    Se = 0 #Se is the leakage area around baffles, assume zero for now
    Sm = d_shell - tube_number * do #maximum free area for flow in cross flow zone
    Sp = d_shell - tube_number * do #minimum free area for flow in cross flow zone

    P = ((pitch - do)/pitch * (do/d_shell))

    Gav = 1/3 (m_dot/Sw + m_dot/Sm + m_dot/Sp)

    R = Re(Gav,do,liquid)

    S = (Sw / (Sw + Se))

    F = (L2 + (L1 - L2)(2*L3/(L1-L2))**0.6)/L1

    Nu = 1.9 * R**0.6 * liquid.Pr**0.3 * P**0.4 * S**2 * F

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

def dP_shell(V,liquid,do,N,tube_layout):

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

def mdot_dP(m_dot,dP_ovr,side,liquid):

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
    

    mdot_from_dP = interp1d(mdot_dP_array[:,1],mdot_dP_array[:,0],fill_value='extrapolate')
    dP_from_mdot = interp1d(mdot_dP_array[:,0],mdot_dP_array[:,1],fill_value='extrapolate')

    dP_new = dP_from_mdot(m_dot/(liquid.rho/1000))
    m_dot_new = mdot_from_dP(dP_ovr)*(liquid.rho/1000)

    rel_e_dP = (dP_new - dP_ovr)/dP_ovr
    rel_e_mdot = (m_dot_new - m_dot)/m_dot

    m_dot = (m_dot + m_dot_new)/2

    return m_dot,rel_e_dP,rel_e_mdot


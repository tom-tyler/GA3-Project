import numpy as np
from iapws._iapws import _Liquid
from hx_classes import HX,water

import hx_functions as hxf
import hydraulic as hyd
import thermal as th


accuracy = 0.001

hx = HX(tube_number = 13,
        baffle_number = 14,
        pitch = 12e-3,
        tube_length = 362e-3,
        shell_length = 450e-3,
        baffle_area = 16e-4,
        tube_layout='t',
        shell_passes=1,
        nozzle_bore=25e-3,
        F=1)

#initial values:
m_h = 0.45 #initial guess for hot mass flow rate
m_c = 0.50 #initial guess for cold mass flow rate

#specified inlet temperatures
T_inh = 53.4
T_inc = 19.7

#need to sort this. curretly issue is _Liquid needs values above 1 bar
dP_measured_c = 0
dP_measured_h = 0

#initial guesses for outlet temperatures
T_outh = 45.3
T_outc = 21

#initial guesses for inlet and outlet pressures
P_inh = 1.01325e5
P_outh = 1.01325e5 - dP_measured_h
P_inc = 1.01325e5
P_outc = 1.01325e5 - dP_measured_c

#initial heat transfer parameters
heat_transfer,eff = 1,1
heat_transfer_new,eff_new = 1,1
Q_counter = 0
rel_e_h1,rel_e_c1 = 1,1

while (abs(rel_e_c1) > accuracy) and (abs(rel_e_c1) > accuracy):

    #creating hot and cold water objects
    h_w = water(T_inh,T_outh,P_inh,P_outh)
    c_w = water(T_inc,T_outc,P_inc,P_outc)

    #HYDRAULIC DESIGN
    hydraulic = hyd.hydraulic_design(m_c,m_h,h_w,c_w,hx,accuracy)

    m_h, m_c = hydraulic['m_h'], hydraulic['m_c']
    dP_hot, dP_cold = hydraulic['dP_hot'], hydraulic['dP_cold']
    V_tube, V_shell = hydraulic['V_tube'], hydraulic['V_shell']
    Ch, Cc = hydraulic['Ch'], hydraulic['Cc']
  
    #THERMAL DESIGN
    thermal = th.thermal_design(Ch,Cc,V_tube,V_shell,hx,h_w,c_w,accuracy,T_inh,T_inc,T_outh,T_outc,hx.F)
    T_outh,T_outc = thermal['T_outh'], thermal['T_outc']
    rel_e_h1,rel_e_c1 = thermal['rel_e_h1'], thermal['rel_e_c1']

    #HEAT TRANSFER

    heat_transfer_h = hxf.Q_h(Ch,T_outh,T_inh)
    heat_transfer_c = hxf.Q_c(Cc,T_outc,T_inc)
    heat_transfer = np.mean([heat_transfer_c,heat_transfer_h])
    eff = hxf.effectiveness(heat_transfer,Cc,Ch,T_inc,T_inh)

    Q_counter += 1

    if Q_counter > 100:
        print('exceeded max iterations for Q')
        break

    #now loop over entire thing again using these 4 values to get new property values to make answer more accurate. when this converges, can find effectiveness and Q
    #need to use lmtd and e-ntu approaches

<<<<<<< HEAD

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
=======
print('Q:',heat_transfer)
print('eff:',eff)
print('T_outh',T_outh,'T_outc',T_outc)
print('m_h',m_h,'m_c',m_c)
print('dP_hot',dP_hot,'dP_cold',dP_cold)
print('mass:',hx.total_mass())
>>>>>>> 56246f44ec8cedee1f9618e5afe5e318b266fc4e

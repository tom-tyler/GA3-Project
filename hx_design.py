import numpy as np
from iapws._iapws import _Liquid
from hx_classes import HX,water

import hx_functions as hxf


accuracy = 0.001

hx = HX(tube_number = 13,
        baffle_number = 14,
        pitch = 12e-3,
        tube_length = 362e-3,
        shell_length = 450e-3,
        baffle_gap = 14e-3,
        baffle_type = 'across_c',
        tube_layout='t',
        shell_passes = 1,
        crossflow_tube_fraction = 1,
        bypass_area_fraction = 0,
        seal_strips = 0,
        crossflow_rows = 4.5,
        tube_bundle_diameter= 56e-3,
        tube_passes = 1,
        year = 2022,
        T_inh = 53.4,
        T_inc = 19.7)

#initial guesses for mass flowrate:
m_h = 0.45 #initial guess for hot mass flow rate
m_c = 0.50 #initial guess for cold mass flow rate

#initial guesses for outlet temperatures
T_outh = 48
T_outc = 22

#initial heat transfer parameters
heat_transfer,eff = 1,1
Q_counter = 0
rel_e_h1,rel_e_c1 = 1,1


while (Q_counter < 10) or ((abs(rel_e_c1) > accuracy) and (abs(rel_e_h1) > accuracy)):

    #creating hot and cold water objects
    h_w = water(hx.T_inh,T_outh)
    c_w = water(hx.T_inc,T_outc)

    #HYDRAULIC DESIGN
    hydraulic = hxf.hydraulic_design(m_c,m_h,h_w,c_w,hx,accuracy)

    m_h, m_c = hydraulic['m_h'], hydraulic['m_c']
    dP_hot, dP_cold = hydraulic['dP_hot'], hydraulic['dP_cold']
    Ch, Cc = hydraulic['Ch'], hydraulic['Cc']
  
    #THERMAL DESIGN
    thermal = hxf.thermal_design(m_h,m_c,h_w,c_w,hx,accuracy,hx.T_inh,hx.T_inc,T_outh,T_outc)
    T_outh,T_outc = thermal['T_outh'], thermal['T_outc']
    rel_e_h1,rel_e_c1 = thermal['rel_e_h1'], thermal['rel_e_c1']

    #HEAT TRANSFER

    heat_transfer_h = hxf.Q_h(Ch,T_outh,hx.T_inh)
    heat_transfer_c = hxf.Q_c(Cc,T_outc,hx.T_inc)
    heat_transfer = np.mean([heat_transfer_c,heat_transfer_h])
    eff = hxf.effectiveness(heat_transfer,Cc,Ch,hx.T_inc,hx.T_inh)

    heat_transfer_ntu = thermal['q_ntu']
    eff_ntu = thermal['eff_ntu']
    U = thermal['U']

    Q_counter += 1
    #print(f'q counter: {Q_counter}')

    if Q_counter > 100:
        print('exceeded max iterations for Q')
        break

    #now loop over entire thing again using these 4 values to get new property values to make answer more accurate. when this converges, can find effectiveness and Q
    #need to use lmtd and e-ntu approaches

print('OUTPUT')
print('--------------------------------------------------------------')
print('HYDRAULIC DESIGN:')
print('\nhot water mass flow rate: ',m_h,'\nhot water pressure drop: ',dP_hot)
print('\ncold water mass flow rate: ',m_c,'\ncold water pressure drop: ',dP_cold)
print('--------------------------------------------------------------')
print('THERMAL DESIGN:')
print('\nLMTD:')
print('T_outh',T_outh,'T_outc',T_outc)
print('Q_LMTD:',heat_transfer)
print('eff_lmtd:',eff)
print('\nNTU:')
print('eff_ntu:',eff_ntu)
print('Q_NTU:',heat_transfer_ntu)
print('--------------------------------------------------------------')
print('OVERALL:')
print('\nmass:',hx.total_mass())
print('U: ',U)
print('--------------------------------------------------------------')
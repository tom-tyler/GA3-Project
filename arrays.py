import numpy as np
import hx_functions as hxf


import numpy as np
from iapws._iapws import _Liquid
from hx_classes import HX,water

import hx_functions as hxf

accuracy = 0.001

hx = HX(13,9,14e-3,350e-3,500e-3,16e-4,tube_layout='s',shell_passes=1)

#initial values:
m_h = 0.6 #initial guess for hot mass flow rate
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

#initial heat transfer parameters
heat_transfer,eff = 1,1
heat_transfer_new,eff_new = 1,1
Q_counter = 0
per_e_Q,per_e_eff = 1,1


#creating hot and cold water objects
h_w = water(T_inh,T_outh,P_inh,P_outh)
c_w = water(T_inc,T_outc,P_inc,P_outc)

m_counter = 0
dP_e = 1
m_e = 1


while (abs(dP_e) > accuracy) and (abs(m_e) > accuracy):

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
      print(dP_tube_ovr)

      # now need iteration routine to get m_h such that dP_tube_ovr matches figure 6 from handout
      print(hxf.mdot_dP(m_h,dP_tube_ovr,'h'))
      m_h, dP_e, m_e = hxf.mdot_dP(m_h,dP_tube_ovr,'h')
      
      m_counter += 1
      print(f'm: {m_h},dP: {dP_tube_ovr}')

      if m_counter > 100:
            print('exceeded max iterations for m,dP')
            break

    

        


mdot_dp_h = [[0.4583,   0.1333],
             [0.4236,   0.1756],
             [0.4010,   0.2024],
             [0.3611,   0.2577],
             [0.3125,   0.3171],
             [0.2639,   0.3633],
             [0.2222,   0.4233],
             [0.1597,   0.4784],
             [0.1181,   0.5330],
             [0.0694,   0.5715]] #[mass flow rate,pressure difference] for hot side

mdot_dp_c = [[0.5833,   0.1113],
             [0.5083,   0.2157],
             [0.4750,   0.2538],
             [0.4250,   0.3168],
             [0.3792,   0.3613],
             [0.3417,   0.4031],
             [0.2958,   0.4511],
             [0.2583,   0.4846],
             [0.2125,   0.5181],
             [0.1708,   0.5573]] #[mass flow rate, pressure difference] for cold side

import hx_functions as hxf
from hx_classes import HX,water
from scipy.optimize import fsolve
import numpy as np

def thermal_design(Ch,Cc,V_tube,V_shell,hx,h_w,c_w,accuracy,T_inh,T_inc,T_outh,T_outc,F):
    
    hi = hxf.hi(V_tube,hx.tube.d_inner,h_w)
    ho = hxf.ho(V_shell,hx.tube.d_outer,c_w,hx.tube_layout)

    U = hxf.U_inside(hi,ho,hx.tube.d_inner,hx.tube.d_outer,hx.tube_length)
    print(U)

    A_con = hx.convection_area

    cmin = min(Cc,Ch)
    cmax = max(Cc,Ch)
    qmax = cmin * (T_inh - T_inc) #maximum possible heat tranfer
    Cr = cmin/cmax #ratio of specific heats 
    NTU = (U * A_con)/cmin
    print(f'counter or co: {hx.co_counter}')
    if  hx.co_counter == 'counter':
        e = (1 - np.exp(-NTU * (1 + Cr)))/(1 + Cr) #equations from wiki, check
    elif hx.co_counter == 'co':
        e = (1 - np.exp(-NTU * (1 - Cr)))/(1 - Cr * np.exp(-NTU * (1 - Cr)))
    else:
        print('Error, heat exchanger must be counter or co flow') 
    #may need something about mixed flow here later
    q_ntu = qmax * e


    #now solve thermal design equations by iteration to get T_outh and T_outc. also find P_outh and P_outc

    T_outc_new,T_outh_new = 1,1
    T_counter = 0
    rel_e_h,rel_e_c = 1,1

    while (abs(rel_e_c) > accuracy) and (abs(rel_e_h) > accuracy):

        T_outc_new = fsolve(lambda T_outc: (T_inc - T_outc) + (1/Cc)*A_con*U*F*hxf.LMTD(T_inc,T_inh,T_outc,T_outh), T_outc)[0]
    
        T_outh_new = fsolve(lambda T_outh: (T_inh - T_outh) - (1/Ch)*A_con*U*F*hxf.LMTD(T_inc,T_inh,T_outc,T_outh), T_outh)[0]
        if T_counter == 0:
            rel_e_c1 = (T_outc_new - T_outc)/T_outc
            rel_e_h1 = (T_outh_new - T_outh)/T_outh
        rel_e_c = (T_outc_new - T_outc)/T_outc
        rel_e_h = (T_outh_new - T_outh)/T_outh
        T_outc = T_outc_new
        T_outh = T_outh_new
        T_counter += 1

        if T_counter > 100:
            print('exceeded max iterations for T')
            break

    thermal = {'T_outh':T_outh,'T_outc':T_outc,'rel_e_c1':rel_e_c1,'rel_e_h1':rel_e_h1,'q_ntu':q_ntu,'eff_ntu':e}
    return thermal


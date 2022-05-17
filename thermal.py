import hx_functions as hxf
from hx_classes import HX,water
from scipy.optimize import fsolve

def thermal_design(Ch,Cc,V_tube,V_shell,hx,h_w,c_w,accuracy,T_inh,T_inc,T_outh,T_outc,F):
    
    hi = hxf.hi(V_tube,hx.tube.d_inner,h_w)
    ho = hxf.ho(V_shell,hx.tube.d_outer,c_w,hx.tube_layout)

    U = hxf.U_inside(hi,ho,hx.tube.d_inner,hx.tube.d_outer,hx.tube_length)
    print(U)

    A_con = hx.convection_area

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

    thermal = {'T_outh':T_outh,'T_outc':T_outc,'rel_e_c1':rel_e_c1,'rel_e_h1':rel_e_h1}
    return thermal


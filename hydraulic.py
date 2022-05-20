import hx_functions as hxf



def hydraulic_design(m_c,m_h,h_w,c_w,hx,accuracy,year):

    m_counter = 0
    dP_e_h = 1
    m_e_h = 1
    dP_e_c = 1
    m_e_c = 1

    #for bore of 25
    K_hmix = 4.9
    K_baffles = 0.53

    #for bore of 20 need to do still
    #K_hmix = 7.0
    #K_baffles = 0.49

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
    
        #pressure drop in a single tube
        dP_tube = hxf.dP_tube(hx.tube_length,hx.tube.d_inner,h_w,V_tube)

        #pressure drop due to separation leaving tubes
        dP_in_plus_out = hxf.dP_inout(h_w,V_tube,sigma)

        #pressure drop due to separation leaving nozzle
        dP_in_plus_out_nozzle = hxf.dP_inout(h_w,V_nozzle_h,sigma_nozzle)

        #head loss in nozzle
        dP_nozzles_h = 2 * hxf.dP_nozzle(V_nozzle_h,h_w)

        #mixing loss
        dP_mixing = K_hmix*h_w.rho*((V_nozzle_h+V_tube)/2)**2

        #overall pressure drop
        dP_tube_ovr = dP_tube + dP_in_plus_out + dP_nozzles_h + dP_in_plus_out_nozzle + dP_mixing

        # now need iteration routine to get m_h such that dP_tube_ovr matches figure 6 from handout
        m_h, dP_e_h, m_e_h = hxf.mdot_dP(m_h,dP_tube_ovr,'h',h_w,year)
      
        #print(f'm_h: {m_h},dP: {dP_tube_ovr}')

        #cold side
        N = hx.tube_number/hx.shell_passes

        dP_shell = hxf.dP_shell(c_w, hx.Ncw, hx.Sb, hx.Sm, hx.Sw, hx.Nc, m_c, hx.baffle_number, hx.b1, hx.b2, hx.b3, hx.b4, hx.tube.d_outer, hx.pitch)
        dP_nozzles_c = 2 * hxf.dP_nozzle(V_nozzle_c,c_w)

        dP_shell_ovr = dP_shell + dP_nozzles_c

        #now need iteration routine to get m_c such that dP_shell_ovr matches figure 6 from handout
        m_c, dP_e_c, m_e_c = hxf.mdot_dP(m_c,dP_shell_ovr,'c',c_w,year)
      
        #print(f'm_c: {m_c},dP: {dP_shell_ovr}')

        m_counter += 1
        if m_counter > 100:
                print('exceeded max iterations for m,dP')
                break

    hydraulic = {'m_h':m_h,'m_c':m_c,'V_tube':V_tube,'V_shell':V_shell,'dP_hot':dP_tube_ovr,'dP_cold':dP_shell_ovr,'Cc':Cc,'Ch':Ch}
    return hydraulic
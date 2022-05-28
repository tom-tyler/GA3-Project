


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


def brute_opt_2():
    #brute optimisation but applying some common sense to reduce time
    baffle_type = 'across_c'
    tube_layout = 't'

    hx_designs = {}
    hx_data = []
    K_hot = 1
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
                                                        name = None
                                                        )
                                    
                                    if heat_exchanger.total_mass() <= 1.1:
                                        print('pass')
                                        hx_designs[f'design {design_no}'] = heat_exchanger
                                        performance = hx_design_basic(heat_exchanger)
                                        
                                        hx_data = hx_data.append(performance, vars(heat_exchanger)[0:13])
    
    hx_data.sort(reverse = True)
    return hx_data[0:3]


def dPw_ideal(m_c, liquid, hx):
    '''Function returns the '''
    dpw_ideal = ((2 + 0.6 * hx.Ncw) * m_c**2) / (2 * 1 * liquid.rho * hx.Sm * hx.Sw)
    return dpw_ideal



def dP_ideal(m_c, liquid, hx):
    #ideal pressure drop
    G = m_c/hx.Sm
    b = hx.b3 / (1 + 0.14 * ((hx.tube.d_outer * G)/liquid.mu)**hx.b4)
    phi = 1
    f_ideal = hx.b1 * (1.33/(hx.pitch/hx.tube.d_outer))**b * ((hx.tube.d_outer * G)/liquid.mu)**hx.b2
    dP_ideal = 
    return dP_ideal



def dP_shell_drop(liquid, m_c, hx, K_cold):
    
    nb = hx.baffle_number
    G = m_c/hx.Sm
    b = hx.b3 / (1 + 0.14 * ((hx.tube.d_outer * G)/liquid.mu)**hx.b4)
    phi = 1
    f_ideal = hx.b1 * (1.33/(hx.pitch/hx.tube.d_outer))**b * ((hx.tube.d_outer * G)/liquid.mu)**hx.b2

    dPi = ((2 + 0.6 * hx.Ncw) * m_c**2) / (2 * 1 * liquid.rho * hx.Sm * hx.Sw)
   
    dPw = (2 * f_ideal * hx.Nc * G**2)/(1 * liquid.rho * phi)
  
    dPs = ((nb - 1)*dPi*hx.Rb + nb*dPw)*hx.Rl*K_cold * hx.shell_passes
    dPe = 2*dPi*(1 + hx.Ncw/hx.Nc)*hx.Rb*hx.Rs 
    dP = dPe + dPs
    return dP

#pressure drop due to separation leaving nozzle
dP_in_plus_out_nozzle = dP_inout(h_w,V_nozzle_h,sigma_nozzle)

def LMTD(T_1in,T_2in,T_1out,T_2out):
    lmtd = ((T_2in - T_1out) - (T_2out - T_1in)) / np.log((T_2in - T_1out) / (T_2out - T_1in))
    return lmtd

# if self.leakage == False:
            # self.Ssb = 0
            # self.Stb = 0
            # self.rs = 1
            # self.rl = 0
            # self.Rl = 1
            # self.Jl = 1 


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

    #NTU1 = NTU/hx.shell_passes
    #e1 = fsolve(lambda e1: NTU1 + (np.log(((2/e1 - (1 + Cr))/c_root - 1)/((2/e1 - (1 + Cr))/c_root + 1)))/c_root, 0.1)[0]

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




    # hx_list['hx_y2017B_p2017'] = HX(tube_number = 14,
    #                     baffle_number = 12,
    #                     pitch = 15e-3,
    #                     tube_length = 225e-3,
    #                     plenum_length_1 = 41e-3,
    #                     plenum_length_2 = 23e-3,
    #                     baffle_gap = 25.6e-3,
    #                     baffle_type = 'across_c',
    #                     tube_layout='t',
    #                     shell_passes = 2,
    #                     tube_bundle_diameter= 60e-3,
    #                     tube_passes = 2,
    #                     baffle_spacing_in = 31e-3,
    #                     baffle_spacing_out = 33e-3,
    #                     design_year = 2017,
    #                     pump_year = 2017,
    #                     T_inh = 52.4,
    #                     T_inc = 21.7,
    #                     name = '2017-B',
    #                     real_data = real_moodle_data['hx_y2017B_p2017'])

    # hx_list['hx_y2017C_p2017'] = HX(tube_number = 14,
    #                     baffle_number = 8,
    #                     pitch = 10e-3,
    #                     tube_length = 192e-3,
    #                     plenum_length_1 = 41e-3,
    #                     plenum_length_2 = 24e-3,
    #                     baffle_gap = 15.75e-3,
    #                     baffle_type = 'across_c',
    #                     tube_layout='t',
    #                     shell_passes = 1,
    #                     tube_bundle_diameter= 52e-3,
    #                     tube_passes = 2,
    #                     baffle_spacing_in = 31.25e-3,
    #                     baffle_spacing_out = 31.25e-3,
    #                     design_year = 2017,
    #                     pump_year = 2017,
    #                     T_inh = 51.5,
    #                     T_inc = 22.2,
    #                     name = '2017-C',
    #                     real_data = real_moodle_data['hx_y2017C_p2017'])


    p_data['hx_y2017B_p2017'] = {'T cold out (C)':26.1,   
                                'T hot out (C)':46.8,
                                'mdot_cold (l/s)':0.454,
                                'mdot_hot (l/s)':0.472,
                                'dP_cold (bar)':0.277,
                                'dP_hot (bar)':0.137,
                                'Q_NTU (kW)':12460,
                                'eff_NTU':0.197,
                                'mass (kg)':1.466
                                }

    p_data['hx_y2017C_p2017'] = {'T cold out (C)':26.1,   
                                'T hot out (C)':46.8,
                                'mdot_cold (l/s)':0.454,
                                'mdot_hot (l/s)':0.472,
                                'dP_cold (bar)':0.277,
                                'dP_hot (bar)':0.137,
                                'Q_NTU (kW)':12460,
                                'eff_NTU':0.197,
                                'mass (kg)':1.466
                                }

    p_data['hx_y2018_p2019'] = {'T cold out (C)':21.0,   
                                'T hot out (C)':45.3,
                                'mdot_cold (l/s)':0.525,
                                'mdot_hot (l/s)':0.490,
                                'dP_cold (bar)':0.356,
                                'dP_hot (bar)':0.163,
                                'Q_NTU (kW)':15820,
                                'eff_NTU':0.198,
                                'mass (kg)':1.466
                                }
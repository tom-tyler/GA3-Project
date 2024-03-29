def ho(do, tube_number, liquid, d_shell, pitch, baffle_area, A_shell, tube_length, thickness, baffle_number, m_dot):
    #relation from https://reader.elsevier.com/reader/sd/pii/0017931063900371?token=67FD7B1EA2E8B9D710BA94CFEC6DA5F06A07655CAAFAFD0A8CC482CCE90F3D4CC6E9878A64AE9F808EAC4F0844192201&originRegion=eu-west-1&originCreation=20220517160904
    L1 = tube_length - 2*thickness #distance between end plates
    L3 = L1/baffle_number - thickness #distance between baffles
    L2 = L1 - 2*L3 - 2*thickness#'length between end baffles' not quite sure how to interpret this

    Sw = A_shell - baffle_area #Sw is the free area for flow to go through in the plane of the baffles
    print(Sw)
    Se = 0 #Se is the leakage area around baffles, assume zero for now
    Sm = d_shell - tube_number * do #maximum free area for flow in cross flow zone
    print(Sm)
    Sp = d_shell - tube_number * do #minimum free area for flow in cross flow zone

    P = ((pitch - do)/pitch * (do/d_shell))

    Gav = 1/3 * (m_dot/Sw + m_dot/Sm + m_dot/Sp)

    R = Re(Gav,do,liquid)

    S = (Sw / (Sw + Se))

    F = (L2 + (L1 - L2)(2*L3/(L1-L2))**0.6)/L1

    Nu = 1.9 * R**0.6 * liquid.Pr**0.3 * P**0.4 * S**2 * F

    ho = Nu*liquid.k/do

    return ho

def ho(V_shell,do,liquid,tube_layout): #handout relation
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
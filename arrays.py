from hx_functions import mdot_dP
from hx_classes import water

h_w = water(60,46)

m = mdot_dP(0.06e5,'h',h_w,2022)
#print(m)


#do not delete what is below!!
                        real_data = {'Name':f'{}-p{()}',
                        'T cold in (C)':,
                        'T cold out (C)':,
                        'T hot in (C)':,
                        'T hot out (C)':,
                        'mdot_cold (l/s)':,
                        'mdot_hot (l/s)':,
                        'dP_cold (bar)':,
                        'dP_hot (bar)':,
                        'Q_LMTD (kW)':,
                        'eff_LMTD':,
                        'Q_NTU (kW)':,
                        'eff_NTU':,
                        'mass (kg)':
                        }

    real_data = {'Name':'JPL-2018-p2022',
                        'T cold in (C)':19.7,
                        'T cold out (C)':26.1,
                        'T hot in (C)':53.4,
                        'T hot out (C)':46.8,
                        'mdot_cold (l/s)':0.454,
                        'mdot_hot (l/s)':0.472
                        'dP_cold (bar)':0.277,
                        'dP_hot (bar)':0.137,
                        'Q_LMTD (kW)':12.5,
                        'eff_LMTD':0.197,
                        'Q_NTU (kW)':12.5,
                        'eff_NTU':0.197,
                        'mass (kg)':1.466}
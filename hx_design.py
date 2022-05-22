import numpy as np
import pandas as pd
import hx_functions as hxf



K_hot = 1.8
K_cold = 1

hx_list = hxf.heat_exchangers()
hx_data = pd.DataFrame()

for hxi in hx_list:
    hx = hx_list[hxi]

    performance = hxf.hx_design(hx,K_hot,K_cold)


    hx_data = hx_data.append(performance, ignore_index = True)

#order columns nicely
hx_data = hx_data[['Name',
               'T cold in (C)',
               'T cold out (C)',
               'T hot in (C)',
               'T hot out (C)',
               'mdot_cold (l/s)',
               'mdot_hot (l/s)',
               'dP_cold (bar)',
               'dP_hot (bar)',
               'Q_LMTD (kW)',
               'eff_LMTD',
               'Q_NTU (kW)',
               'eff_NTU',
               'mass (kg)'
                ]]
with pd.option_context('display.max_rows', None, 'display.max_columns', None,"display.precision", 3):  # more options can be specified also
    print(hx_data)
import matplotlib.pyplot as plt
from hx_functions import hx_design, fit_data
from hx_classes import HX

performance = []
tn = []

k_array = [5.09582324e-13, 1.06322974e-02, 1.60725045e+01, 7.70041192e-14,
9.74222768e-01, 7.39759976e+00, 6.78054364e+01, 7.84344345e-01,
3.25600821e-01, 8.46334898e-01]

for n in range(1,200):
    hx =  HX(tube_number = 20,
             baffle_number = 6,
             pitch = 10e-3,
             tube_length = 172e-3,
             plenum_length_1 = 58e-3,
             plenum_length_2 = 23e-3,
             baffle_gap = n/10000,
             baffle_type = 'across_c',
             tube_layout='t',
             shell_passes = 2,
             tube_bundle_diameter= 54e-3,
             tube_passes = 4,
             baffle_spacing_in = 41.5e-3,
             baffle_spacing_out = 21.5e-3,
             design_year = 2022,
             pump_year = 2022,
             T_inh = 60,
             T_inc = 20,
             name = '2022_F',
             real_data = None)

    tn.append(n)
    performance.append(hx_design(hx,k_array)['Q_NTU (kW)'])

plt.plot(tn,performance)
plt.show()


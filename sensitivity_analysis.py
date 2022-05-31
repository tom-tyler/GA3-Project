import matplotlib.pyplot as plt
from hx_functions import hx_design, fit_data
from hx_classes import HX

performance = []
tn = []

for tube_number in range(15,40):
    hx =  HX(tube_number = tube_number,
             baffle_number = 6,
             pitch = 10e-3,
             tube_length = 172e-3,
             plenum_length_1 = 58e-3,
             plenum_length_2 = 23e-3,
             baffle_gap = 20e-3,
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

    k_array = fit_data()
    tn.append(tube_number)
    performance.append(hx_design(hx)['Q_NTU (kW)'])

plt.plot(tn,performance)
plt.show()


import matplotlib.pyplot as plt
from hx_functions import hx_design, fit_data
from hx_classes import HX

performance = []
tn = []

k_array = [0.75,       0.75,       2.,         0.5,        1.52671187, 1.25,
1.5,        1.5,        0.60125398, 0.9503981] #once this is done once just copy the array to save time

for n in range(1,20):
    hx =  HX(tube_number=20,
             baffle_number=6,
             pitch=10e-3,
             tube_length=172e-3,
             plenum_length_1=58e-3,
             plenum_length_2=23e-3,
             baffle_gap=20.00e-3,
             baffle_type='across_c',
             tube_layout='t',
             shell_passes=2,
             tube_bundle_diameter=55e-3,
             tube_passes=4,
             design_year=2022,
             pump_year=2022,
             T_inh=60,
             T_inc=20,
             baffle_spacing_in=41.5e-3,
             baffle_spacing_out=21.5e-3,
             name='2022_F',
             real_data=None
             )

    tn.append(n)
    performance.append(hx_design(hx,k_array)['Q_NTU (kW)'])

# let us make a simple graph
fig = plt.figure(figsize=[7,5])
ax = plt.subplot(111)
plt.plot(tn,performance)
# l = ax.fill_between(tn,performance)

# set the basic properties
ax.set_xlabel('Baffle Gap (mm)')
ax.set_ylabel('Q (W)')
# ax.set_title('Title')

# set the limits
ax.set_xlim(0, 20)
ax.set_ylim(10, 20)

# set the grid on
ax.grid('on')

plt.show()

# plt.rcParams["figure.figsize"] = (5,5)
# plt.plot(tn,performance)
# plt.ylabel('Baffle Gap (mm)')
# plt.xlabel('Q (W)')
# # plt.title('Worked Example Cold Hydraulic Design')
# # plt.legend()
# plt.show()


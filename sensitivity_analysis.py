import matplotlib.pyplot as plt
from hx_functions import hx_design
from hx_classes import HX
import numpy as np

performance1 = []
tn1 = []
performance2 = []
tn2 = []
performance3 = []
tn3 = []

k_array = [0.75,       0.75,       2.,         0.5,        1.52671187, 1.25,
1.5,        1.5,        0.60125398, 0.9503981] #once this is done once just copy the array to save time

for n in np.arange(5,25,1):
    hx =  HX(tube_number=n,
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

    tn1.append(n)
    performance1.append(hx_design(hx,k_array)['Q_NTU (kW)'])

for n in np.arange(2,20,1):
    hx =  HX(tube_number=20,
             baffle_number=n,
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

    tn2.append(n)
    performance2.append(hx_design(hx,k_array)['Q_NTU (kW)'])

for n in np.linspace(5e-3, 20e-3):
    hx =  HX(tube_number=20,
             baffle_number=6,
             pitch=10e-3,
             tube_length=172e-3,
             plenum_length_1=58e-3,
             plenum_length_2=23e-3,
             baffle_gap=n,
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

    tn3.append(n)
    performance3.append(hx_design(hx,k_array)['Q_NTU (kW)'])

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 12}
plt.rc('font', **font)
fig, axs = plt.subplots(1,3,figsize=(16,5), constrained_layout=True)
fig.suptitle('Sensitivity Analysis')
axs[0].plot(tn1,performance1,marker='o', c='turquoise')
axs[0].set_xlabel('Tube Number')
axs[0].set_ylabel('$Q_{corr}$ (kW)')
axs[0].set_xlim([0, 30])

axs[1].plot(tn2, performance2, marker='o', c='firebrick')
axs[1].set_xlabel('Baffle Number')
axs[1].set_ylabel('$Q_{corr}$ (kW)')
axs[1].set_xlim([0,20])

axs[2].plot(tn3, performance3, c='violet')
axs[2].set_xlabel('Baffle Gap (m)')
axs[2].set_ylabel('$Q_{corr}$ (kW)')
axs[2].set_xlim([0,25e-3])

plt.show()
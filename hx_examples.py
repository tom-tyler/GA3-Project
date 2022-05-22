

#2017 A1
hx = HX(tube_number = 8,
        baffle_number = 12,
        pitch = 10e-3,
        tube_length = 261e-3,
        shell_length = 334e-3,
        baffle_gap = 16e-3,
        baffle_type = 'across',
        tube_layout='t',
        shell_passes=1,
        tube_passes=1,
        nozzle_bore=19e-3)

#jpl 2018
hx = HX(tube_number = 13,
        baffle_number = 14,
        pitch = 12e-3,
        tube_length = 362e-3,
        shell_length = 450e-3,
        baffle_gap = 14e-3,
        baffle_type = 'across_c',
        tube_layout='t',
        shell_passes=1,
        nozzle_bore=25e-3)


hx = HX(tube_number = 10,
           baffle_number = 10,
           pitch = 10/1000,
           tube_length = 200e-3,
           shell_length = 300e-3,
           baffle_gap = 10e-3,
           baffle_type = 'across_c',
           tube_layout='t',
           shell_passes=1,
           nozzle_bore=20e-3,
           crossflow_tube_fraction = 1,
           bypass_area_fraction = 0,
           seal_strips = 0,
           crossflow_rows = 4.5,
           tube_bundle_diameter= 56e-3,
           tube_passes = 1)

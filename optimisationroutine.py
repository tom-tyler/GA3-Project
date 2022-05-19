# first, brute force optimisation
import hx_functions as hxf
from hx_classes import HX
n = 0



for tube_length in range(0,269):#incorporates length constraint
    for tube_number in range(1,30):
        if tube_length * tube_number <= 3500:
            #run design code
            print('solution')
            for baffle_number in range(0,50):
                for shell_passes in range(1,4):
                    hx = HX(tube_number = tube_number,
                            baffle_number = baffle_number,
                            pitch = 12e-3,
                            tube_length = tube_length,
                            shell_length = 450e-3,
                            baffle_gap  = 13.9e-3,
                            baffle_type = 'across',
                            tube_layout='t',
                            shell_passes = shell_passes,
                            tube_passes=1,
                            nozzle_bore=25e-3)
                    if HX.total_mass(HX) <= 1.1:
                        print('fair solution')
                        n+=1

print(n)

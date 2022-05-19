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
                    if HX.total_mass(HX) <= 1.1:
                        print('fair solution')
                        n+=1

print(n)

from geneal.genetic_algorithms import ContinuousGenAlgSolver
from geneal.applications.fitness_functions.continuous import fitness_functions_continuous
from hx_classes import HX
import hx_functions as hxf
import hx_design as hxd
import numpy as np
import matplotlib.pyplot as plt

solver = ContinuousGenAlgSolver(n_genes=4, # number of variables defining the problem
                                fitness_function=fitness_functions_continuous(3), # fitness function to be maximized
                                pop_size=10,
                                max_gen=200,
                                mutation_rate=0.1,
                                selection_rate=0.6,
                                selection_strategy="roulette_wheel",
                                problem_type=float, # Defines the possible values as float numbers
                                variables_limits=(-10, 10) # Defines the limits of all variables between -10 and 10. 
                                                           # Alternatively one can pass an array of tuples defining the limits
                                                           # for each variable: [(-10, 10), (0, 5), (0, 5), (-20, 20)]
                                )

#solver.solve()


def brute_opt(hx,step):
    #brute force optimisation with a few checks to eliminate cases as early as possible

    baffle_type = 'across_c'
    tube_layout = 't'

    designs = []

    for tube_passes in range(1,4):
        for tube_number in range(12,15):
            #print(tube_number) #just gives an idea of progress
            #pitch_max = round(hx.shell.d_inner / tube_number)
            if tube_passes%2 == 0:
                l_min = 1
            else:
                l_min = 41
            for plenum_length_1 in range(41,100, step):
                for plenum_length_2 in range(l_min, 100,step):
                    for tube_length in range(150,250,step):
                        if tube_number * tube_length <= 3500:
                            for shell_passes in range(1,4):
                                if shell_passes%2 == 0:
                                    bso_min = 1
                                else:
                                    bso_min = 41
                                for baffle_spacing_in in range(41,100,step):
                                    for baffle_spacing_out in range(bso_min, 100, step):
                                        for baffle_number in range(10,13):
                                            for baffle_gap in range(10,60,10):
                                                #can we do mass constraint here
                                                    for pitch in range(10,20,10):
                                                        for crossflow_rows in range(4,5):
                                                            hx = HX(tube_number = tube_number,
                                                                    baffle_number = baffle_number,
                                                                    pitch = pitch/1000,
                                                                    tube_length = tube_length/1000,
                                                                    plenum_length_1 = plenum_length_1/1000,
                                                                    plenum_length_2 = plenum_length_2/1000,
                                                                    baffle_gap = baffle_gap/1000,
                                                                    baffle_type = baffle_type,
                                                                    tube_layout = tube_layout,
                                                                    shell_passes = shell_passes,
                                                                    tube_bundle_diameter  = (crossflow_rows * pitch + 8e-3)/1000,
                                                                    tube_passes = tube_passes,
                                                                    baffle_spacing_in = baffle_spacing_in/1000,
                                                                    baffle_spacing_out = baffle_spacing_out/1000,
                                                                    design_year = 2022,
                                                                    pump_year = 2022,
                                                                    T_inh = 53.4,
                                                                    T_inc = 19.2,
                                                                    leakage = True,
                                                                    name = None,
                                                                    co_counter='counter',
                                                                    approximate_glue_mass=0
                                                                    )
                                                            if hxf.total_mass(hx) <= 1.1:
                                                                q = hxd.hx_design(hx)
                                                                if np.isnan(q) == False:
                                                                    attributes = vars(hx)
                                                                    designs.append((q,list(attributes.items())[0:13]))


    designs.sort(reverse = True)
    return designs[0:2], tube_length


            
def optimiser(hx):
    best_array = []
    step_array = []
    tube_length_array = []
    step = 50

    optimised = brute_opt(hx,step)
    best_design = optimised[0]
    q_best = best_design[0]
    tube_length = optimised[1]
    best_array.append(q_best[0])
    tube_length_array.append(tube_length*20)
    step_array.append(step)

    #plt.plot(step_array, best_array)
    #plt.plot(step_array, tube_length_array)
    plt.show()

    return optimised


print(optimiser(HX))



#brute_opt(HX)


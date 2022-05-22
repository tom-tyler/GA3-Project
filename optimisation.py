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

    for tube_number in range(12,15):
        #print(tube_number) #just gives an idea of progress
        #pitch_max = round(hx.shell.d_inner / tube_number)
        for shell_length in range(300,350, step):
            #print(shell_length)
            for tube_length in range(150,250,step):
                #print(tube_length)
                if shell_length - tube_length >= 10: #check this constraint, need to be able to fit nozzles on
                    if tube_number * tube_length <= 3500:
                        for shell_passes in range(1,4):
                            for baffle_number in range(10,13):
                                for baffle_gap in range(10,60,10):
                                    for seal_strips in range(0,1):
                                        for pitch in range(10,20,10):
                                            for tube_passes in range(1,4):
                                                for crossflow_rows in range(4,5):
                                                    hx = HX(tube_number = tube_number,
                                                            baffle_number = baffle_number,
                                                            pitch,
                                                            tube_length,
                                                            plenum_length_1,
                                                            plenum_length_2,
                                                            baffle_gap,
                                                            baffle_type,
                                                            tube_layout,
                                                            shell_passes,
                                                            tube_bundle_diameter,
                                                            tube_passes,
                                                            baffle_spacing_in,
                                                            baffle_spacing_out,
                                                            design_year,
                                                            pump_year,
                                                            T_inh,
                                                            T_inc,
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
    for step in range(100,5,-5):
        print(step)
        optimised = brute_opt(hx,step)
        best_design = optimised[0]
        q_best = best_design[0]
        tube_length = optimised[1]
        best_array.append(q_best[0])
        tube_length_array.append(tube_length*20)
        step_array.append(step)

    plt.plot(step_array, best_array)
    plt.plot(step_array, tube_length_array)
    plt.show()

    return optimised


print(optimiser(HX))



#brute_opt(HX)

